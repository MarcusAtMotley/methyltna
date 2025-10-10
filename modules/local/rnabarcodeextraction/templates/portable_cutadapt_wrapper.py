#!/usr/bin/env -S uv run --script
#./
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "cutadapt",
#   "pyyaml",
#   "tqdm",
#   "biopython",
#   "loguru",
# ]
# ///
"""
portable_cutadapt_wrapper.py
Marcus Viscardi, 8/19/2025

General goal of this script is to build a simple wrapper of cutadapt and leverage
uv so that it can be run on different devices with ease!

We'll want defaults for the parser based on our basic NNSR and mNNSR barcodes,
and then allow for some flexibility in the command line arguments.

Let's use a yml file to store the defaults, so that we can easily change them
as we use different barcode options!
"""
from sys import executable
from pathlib import Path
from typing import Tuple
from json import load
from subprocess import run

import tqdm.auto as tqdm
from Bio import SeqIO

from tempfile import NamedTemporaryFile

from loguru import logger
import sys

VALID_LOG_LEVELS = ["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"]

class PairedEndCutAdaptRunner:
    def __init__(self,
                 fastq_r1: Path,
                 fastq_r2: Path,
                 output_dir: Path,
                 output_name_override: str = None,
                 adapter_name: str = None,
                 adapter_sequence: str = None,
                 error_rate: float = 0.0,
                 reverse_complement_search: bool = False,
                 min_overlap: int = 6,
                 action: str = "trim",
                 flip_reads: bool = False,
                 subsample: int = -1,
                 adapter_name2: str = None,
                 adapter_sequence2: str = None,
                 times: int = 1,
                 **kwargs
                 ):

        self.flip_reads = flip_reads
        if flip_reads:
            self.read_1 = fastq_r2
            self.read_2 = fastq_r1
        else:
            self.read_1 = fastq_r1
            self.read_2 = fastq_r2
        self.adapter_name = adapter_name if adapter_name else "Adapter"
        self.adapter_sequence = adapter_sequence
        assert adapter_sequence is not None, "Adapter sequence must be provided."
        self.adapter_sequence2 = adapter_sequence2
        self.adapter_name2 = adapter_name2 if adapter_name2 else "Adapter2"
        self.times = times
        if adapter_sequence2 is not None and self.times <= 1:
            self.times = 2
            logger.warning("adapter_sequence2 provided but times param <= 1. Setting times param to 2.")
        if self.adapter_sequence2 is not None:
            self.dual_adapters = True
            assert adapter_name != adapter_name2
        else:
            self.dual_adapters = False
        self.error_rate = error_rate
        self.reverse_complement_search = reverse_complement_search
        self.min_overlap = min_overlap
        self.action = action
        assert self.action in ("trim", "retain", "none", "mask", "crop", "lowercase")
        self.output_dir = output_dir
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)
        if output_dir is None:
            logger.info(f"Output directory not provided. Going to try to make a subdirectory called 'cutadapt' next to the R1 read.")
            self.output_dir = self.read_1.parent / "cutadapt"
        if not self.output_dir.exists():
            if not kwargs.get("do_not_make_output_dir", False):
                self.output_dir.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError(f"Output directory {self.output_dir} does not exist and "
                                        "do_not_make_output_dir is set to True. Please create the directory manually.")
        if not self.output_dir.exists():
            # Weird situtaion where the output directory does not exist, but we are not allowed to create it.
            # Roll back to the parent directory of the first read.
            logger.warning(f"Output directory {self.output_dir} does not exist and could not be created(?) "
                  f"Using parent directory of R1 read instead.")
            self.output_dir = self.read_1.parent
        assert self.output_dir.exists(), f"Output directory {self.output_dir} does not exist. Please create it manually."

        self.output_name_override = output_name_override
        if not self.output_name_override:
            logger.debug("No output name override, using default name.")
            self.output_file_r1 = self.output_dir / self.read_1.name.replace(".fastq", ".cutadapt.fastq")
            self.output_file_r2 = self.output_dir / self.read_2.name.replace(".fastq", ".cutadapt.fastq")
            self.json_report = self.output_dir / f"{self.read_1.name}.cutadapt.json"
            self.txt_report = self.output_dir / f"{self.read_1.name}.cutadapt.txt"
        else:
            logger.debug(f"Output name override: {self.output_name_override}")
            self.output_file_r1 = self.output_dir / f"{self.output_name_override}_R1.cutadapt.fastq"
            self.output_file_r2 = self.output_dir / f"{self.output_name_override}_R2.cutadapt.fastq"
            self.json_report = self.output_dir / f"{self.output_name_override}.cutadapt.json"
            self.txt_report = self.output_dir / f"{self.output_name_override}.cutadapt.txt"

        # This needs to be last so the output files are set correctly!
        if subsample > 0:
            max_pairs = subsample
            tmp_r1 = NamedTemporaryFile(delete=False, suffix=".fastq")
            tmp_r2 = NamedTemporaryFile(delete=False, suffix=".fastq")
            logger.info(f"Subsampling {max_pairs} pairs from R1/R2 reads.")
            with open(self.read_1, "r") as r1, open(self.read_2, "r") as r2, \
                 open(tmp_r1.name, "w") as out_r1, open(tmp_r2.name, "w") as out_r2:
                for i in range(max_pairs):
                    for _ in range(4):
                        out_r1.write(r1.readline())
                    for _ in range(4):
                        out_r2.write(r2.readline())
            self.read_1 = Path(tmp_r1.name)
            self.read_2 = Path(tmp_r2.name)
        self.double_tagged_reads = -1
        self.single_tagged_reads = -1
        self.untagged_reads = -1
        self.frac_tagged_reads = -1


    def run(self, threads=1, extra_params=None, quiet=False, **kwargs) -> Tuple[str, str]:
        cutadapt_call = [
            executable, "-m", "cutadapt",
            f"--action={self.action}",
            "-j", str(threads),
        ]
        # For read1
        cutadapt_call.extend(["-g", f"{self.adapter_name}={self.adapter_sequence}"])
        # For read2
        cutadapt_call.extend(["-G", f"{self.adapter_name}={self.adapter_sequence}"])
        if self.adapter_sequence2 is not None:
            cutadapt_call.extend(["-g", f"{self.adapter_name2}={self.adapter_sequence2}"])
            cutadapt_call.extend(["-G", f"{self.adapter_name2}={self.adapter_sequence2}"])
        add_to_cutadapt_call = [
            "--overlap", str(self.min_overlap),
            "--error-rate", str(self.error_rate),
            "--rename", "{id} Barcode={r1.adapter_name}:{r2.adapter_name} Seq={r1.match_sequence}:{r2.match_sequence} {comment}",
            "--minimum-length", "1",  # Discard reads with zero length after trimming
            str(self.read_1), str(self.read_2),
            "-o", str(self.output_file_r1), "-p", str(self.output_file_r2),
            f"--json={self.json_report}",
            # f"--wildcard-file={self.wildcard_report}",
        ]
        cutadapt_call.extend(add_to_cutadapt_call)
        if self.times != 1:
            cutadapt_call.extend(["--times", str(self.times)])
        if self.reverse_complement_search:
            cutadapt_call.append("--rc")
        if isinstance(extra_params, list) or isinstance(extra_params, tuple):
            cutadapt_call.extend(extra_params)
        elif isinstance(extra_params, str):
            cutadapt_call.append(extra_params)
        if quiet:
            cutadapt_call.append("--quiet")
        # print("Running cutadapt with the following command:")
        # print(" ".join(cutadapt_call))
        logger.debug(f"Running cutadapt with the following command: {' '.join(cutadapt_call)}")
        result = run(cutadapt_call, capture_output=True, text=True)
        logger.debug(f"Cutadapt completed with return code {result.returncode}. ")
        if result.returncode != 0:
            raise RuntimeError(f"Cutadapt failed: {result.stderr}")

        # Clean up sample names in cutadapt report for MultiQC compatibility
        # MultiQC concatenates paired-end input filenames, creating doubled-up names like:
        # "Sample_1_val_1.fq Sample_2_val_2.fq" -> clean to just "Sample.fq Sample.fq"
        import re
        cleaned_stdout = result.stdout
        # Pattern matches: "Sample_1_val_1.fq Sample_2_val_2.fq" or similar
        # Replaces with just the base sample name (first capture group)
        cleaned_stdout = re.sub(r'(\S+?)_1(?:_val_1)?\.f(?:ast)?q\s+\1_2(?:_val_2)?\.f(?:ast)?q',
                               r'\1.fq \1.fq',  # Keep .fq extension for both files
                               cleaned_stdout)

        with open(self.txt_report, "w") as f:
            f.write(cleaned_stdout)
        return result.stdout, result.stderr

    def read_json_report(self) -> dict:
        if not self.json_report.exists():
            raise FileNotFoundError(f"JSON report {self.json_report} does not exist.")
        with open(self.json_report, 'r') as f:
            return load(f)

    def print_quick_report(self):
        # TODO: check if this is broken by dual adapters!!
        if not self.json_report.exists():
            raise ValueError("JSON report not set. Run the cutadapt command first. Or use read_json_report() to load it.")
        report = self.read_json_report()
        read_counts = report.get('read_counts', {})
        total_out_reads = read_counts.get('output', 0)
        total_found_adapters_r1 = read_counts.get('read1_with_adapter', 0)
        total_found_adapters_r2 = read_counts.get('read2_with_adapter', 0)
        # Lets turn this whole thing into a single string for easy logging
        quick_report = (
            f"Quick Report:\n"
            f"\tInput files: {self.read_1}, {self.read_2}\n"
            f"\tOutput files: {self.output_file_r1}, {self.output_file_r2}\n"
            f"\tAdapter sequence: {self.adapter_sequence}\n"
            f"\tTotal output reads:       {total_out_reads:>9} reads\n"
            f"\tTotal reads with adapter R1: {total_found_adapters_r1:>9} "
            f"reads ({total_found_adapters_r1 / total_out_reads:.2%})\n"
            f"\tTotal reads with adapter R2: {total_found_adapters_r2:>9} "
            f"reads ({total_found_adapters_r2 / total_out_reads:.2%})\n"
            f"\tTotal adapters found: {total_found_adapters_r1 + total_found_adapters_r2:>9} reads "
            f"({(total_found_adapters_r1 + total_found_adapters_r2) / total_out_reads:.2%})\n"
            f"\tJSON report: {self.json_report}\n"
            f"\tText report: {self.txt_report}\n"
        )
        logger.info(quick_report)

    def provide_mini_results_dict(self):
        """
        Provides a mini results dictionary with basic information about the run.
        """
        if not self.json_report.exists():
            raise ValueError("JSON report not set. Run the cutadapt command first. Or use read_json_report() to load it.")
        report = self.read_json_report()
        read_counts = report.get('read_counts', {})
        total_out_reads = read_counts.get('output', 0)
        total_found_adapters_r1 = read_counts.get('read1_with_adapter', 0)
        total_found_adapters_r2 = read_counts.get('read2_with_adapter', 0)

        return {
            "output_file_r1": str(self.output_file_r1),
            "output_file_r2": str(self.output_file_r2),
            "adapter_sequence": self.adapter_sequence,
            "total_output_reads": total_out_reads,
            "total_found_adapters_r1": total_found_adapters_r1,
            "total_found_adapters_r2": total_found_adapters_r2,
            "percentage_found_r1": (total_found_adapters_r1 / total_out_reads * 100) if total_out_reads > 0 else 0,
            "percentage_found_r2": (total_found_adapters_r2 / total_out_reads * 100) if total_out_reads > 0 else 0,
            "overall_percentage_barcoded": self.frac_tagged_reads,
            "number_tagged": self.single_tagged_reads + self.double_tagged_reads,
            "number_double_tagged": self.double_tagged_reads,
            "number_single_tagged": self.single_tagged_reads,
            "number_untagged": self.untagged_reads,
            "whole_output_json": report,
        }

    def update_json_report(self):
        # We want to load the existing json report, update it with the new information, and then write it back out.
        if not self.json_report.exists():
            raise ValueError("JSON report not set. Run the cutadapt command first. Or use read_json_report() to load it.")
        report = self.read_json_report()
        add_on_info = {
            "adapter_name": self.adapter_name,
            "adapter_sequence": self.adapter_sequence,
            "output_file_r1": str(self.output_file_r1),
            "output_file_r2": str(self.output_file_r2),
            "number_double_tagged": self.double_tagged_reads,
            "number_single_tagged": self.single_tagged_reads,
            "number_untagged": self.untagged_reads,
            "overall_percentage_barcoded": self.frac_tagged_reads,
            "number_tagged": self.single_tagged_reads + self.double_tagged_reads,
        }
        report.update({"additional_info": add_on_info})
        with open(self.json_report, 'w') as f:
            import json
            json.dump(report, f, indent=2)
        logger.success(f"JSON report updated successfully.")

    def split_tagged_and_untagged_reads(self,
                                        name_for_hits: str = "tagged",
                                        name_for_misses: str = "untagged",
                                        ) -> None:
        """
        Splits the reads into tagged and untagged based on the presence of adapter sequence.
        """
        import gzip
        from mimetypes import guess_type
        from functools import partial

        encoding = guess_type(self.output_file_r1)[1]  # uses file extension
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

        # Create new filename stems with barcoded info embedded in the main filename
        # Instead of: sample.cutadapt_barcoded.fastq (which causes methylSNP collisions)
        # Create: sample_barcoded.cutadapt.fastq (so basename includes the barcode designation)

        # Handle both override and non-override cases properly
        if self.output_name_override:
            # For override case: "override_R1.cutadapt.fastq" -> "override_barcoded_R1.cutadapt.fastq"
            base_name = self.output_name_override
            tagged_fastq_r1 = self.output_file_r1.parent / f"{base_name}_{name_for_hits}_R1.cutadapt.fastq"
            tagged_fastq_r2 = self.output_file_r2.parent / f"{base_name}_{name_for_hits}_R2.cutadapt.fastq"
            untagged_fastq_r1 = self.output_file_r1.parent / f"{base_name}_{name_for_misses}_R1.cutadapt.fastq"
            untagged_fastq_r2 = self.output_file_r2.parent / f"{base_name}_{name_for_misses}_R2.cutadapt.fastq"
        else:
            # For non-override case: "sample_R1_001.cutadapt.fastq" -> "sample_barcoded_R1_001.cutadapt.fastq"
            # Extract the base name before .cutadapt
            r1_name_no_ext = self.read_1.name.replace(".fastq", "").replace(".fq", "")
            r2_name_no_ext = self.read_2.name.replace(".fastq", "").replace(".fq", "")
            tagged_fastq_r1 = self.output_file_r1.parent / f"{r1_name_no_ext}_{name_for_hits}.cutadapt.fastq"
            tagged_fastq_r2 = self.output_file_r2.parent / f"{r2_name_no_ext}_{name_for_hits}.cutadapt.fastq"
            untagged_fastq_r1 = self.output_file_r1.parent / f"{r1_name_no_ext}_{name_for_misses}.cutadapt.fastq"
            untagged_fastq_r2 = self.output_file_r2.parent / f"{r2_name_no_ext}_{name_for_misses}.cutadapt.fastq"
        with _open(self.output_file_r1) as handle_r1:
            total_reads = sum(1 for _ in SeqIO.parse(handle_r1, "fastq"))
        with open(tagged_fastq_r1, 'w') as tagged_r1,\
             open(tagged_fastq_r2, 'w') as tagged_r2,\
             open(untagged_fastq_r1, 'w') as untagged_r1,\
             open(untagged_fastq_r2, 'w') as untagged_r2, \
             _open(self.output_file_r1) as handle_r1, \
             _open(self.output_file_r2) as handle_r2:
            iterator = tqdm.tqdm(enumerate(zip(SeqIO.parse(handle_r1, "fastq"),
                                               SeqIO.parse(handle_r2, "fastq"))),
                                 total=total_reads)
            adapter_names = {self.adapter_name}
            if self.dual_adapters:
                adapter_names.add(self.adapter_name2)
            single_tags, double_tags, untags = 0, 0, 0
            for i, (record_r1, record_r2) in iterator:
                comments_r1 = record_r1.description.split()
                comments_r2 = record_r2.description.split()
                assert comments_r1[1:2] == comments_r2[1:2], \
                    (f"Comments do not match between R1 and R2: {comments_r1} vs {comments_r2}. "
                     f"Could read order have been changed?")
                barcode_1, barcode_2 = comments_r1[1].split("=")[1].split(":")
                if barcode_1 in adapter_names and barcode_2 in adapter_names:
                    # This is a double tagged read (interesting? wouldn't expect this if you only have one adapter!)
                    SeqIO.write(record_r1, tagged_r1, "fastq")
                    SeqIO.write(record_r2, tagged_r2, "fastq")
                    double_tags += 1
                elif barcode_1 in adapter_names or barcode_2 in adapter_names:
                    # This is a tagged read (one of the reads has the adapter)
                    SeqIO.write(record_r1, tagged_r1, "fastq")
                    SeqIO.write(record_r2, tagged_r2, "fastq")
                    single_tags += 1
                else:
                    # This is an untagged read
                    SeqIO.write(record_r1, untagged_r1, "fastq")
                    SeqIO.write(record_r2, untagged_r2, "fastq")
                    untags += 1
                if (i % 100 == 0 and i > 0) or i == len(iterator)-1:
                    iterator.set_description(f"Splitting reads: "
                                             f"{name_for_hits.title()}: {single_tags} ({single_tags / i:.1%}) "
                                             f"Double {name_for_hits.title()}: {double_tags} ({double_tags / i:.1%}) "
                                             f"{name_for_misses.title()}: {untags} ({untags / i:.1%})")
        logger.debug(f"Tagged reads written to {tagged_fastq_r1} and {tagged_fastq_r2}.")
        logger.debug(f"Untagged reads written to {untagged_fastq_r1} and {untagged_fastq_r2}.")
        logger.success(f"{name_for_hits.title()}: {single_tags} ({single_tags / i:.1%}) "
                    f"Double {name_for_hits.title()}: {double_tags} ({double_tags / i:.1%}) "
                    f"{name_for_misses.title()}: {untags} ({untags / i:.1%})")
        # Now we should delete the original output files, since we have split them
        try:
            self.output_file_r1.unlink(missing_ok=True)
            self.output_file_r2.unlink(missing_ok=True)
        except Exception as e:
            logger.info(f"Could not delete original output files: {e}")
            logger.info("You may want to delete them manually if you don't need them anymore.")
        self.double_tagged_reads = double_tags
        self.single_tagged_reads = single_tags
        self.untagged_reads = untags
        self.frac_tagged_reads = (single_tags + double_tags) / (single_tags + double_tags + untags) if (single_tags + double_tags + untags) > 0 else 0

class SingleEndCutAdaptRunner:
    def __init__(self,
                 fastq: Path,
                 output_dir: Path,
                 output_name_override: str = None,
                 adapter_name: str = None,
                 adapter_sequence: str = None,
                 error_rate: float = 0.0,
                 reverse_complement_search: bool = False,
                 min_overlap: int = 6,
                 action: str = "trim",
                 subsample: int = -1,
                 adapter_name2: str = None,
                 adapter_sequence2: str = None,
                 times: int = 1,
                 **kwargs
                 ):

        self.read = fastq
        self.adapter_name = adapter_name if adapter_name else "Adapter"
        self.adapter_sequence = adapter_sequence
        assert adapter_sequence is not None, "Adapter sequence must be provided."
        self.adapter_sequence2 = adapter_sequence2
        self.adapter_name2 = adapter_name2 if adapter_name2 else "Adapter2"
        self.times = times
        if adapter_sequence2 is not None and self.times <= 1:
            self.times = 2
            logger.warning("adapter_sequence2 provided but times param <= 1. Setting times param to 2.")
        if self.adapter_sequence2 is not None:
            self.dual_adapters = True
            assert adapter_name != adapter_name2
        else:
            self.dual_adapters = False
        self.error_rate = error_rate
        self.reverse_complement_search = reverse_complement_search
        self.min_overlap = min_overlap
        self.action = action
        assert self.action in ("trim", "retain", "none", "mask", "crop", "lowercase")
        self.output_dir = output_dir
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)
        if output_dir is None:
            logger.info(f"Output directory not provided. Going to try to make a subdirectory called 'cutadapt' next to the read.")
            self.output_dir = self.read.parent / "cutadapt"
        if not self.output_dir.exists():
            if not kwargs.get("do_not_make_output_dir", False):
                self.output_dir.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError(f"Output directory {self.output_dir} does not exist and "
                                        "do_not_make_output_dir is set to True. Please create the directory manually.")
        if not self.output_dir.exists():
            # Weird situation where the output directory does not exist, but we are not allowed to create it.
            # Roll back to the parent directory of the read.
            logger.warning(f"Output directory {self.output_dir} does not exist and could not be created(?) "
                  f"Using parent directory of read instead.")
            self.output_dir = self.read.parent
        assert self.output_dir.exists(), f"Output directory {self.output_dir} does not exist. Please create it manually."

        self.output_name_override = output_name_override
        if not self.output_name_override:
            logger.debug("No output name override, using default name.")
            self.output_file = self.output_dir / self.read.name.replace(".fastq", ".cutadapt.fastq")
            self.json_report = self.output_dir / f"{self.read.name}.cutadapt.json"
            self.txt_report = self.output_dir / f"{self.read.name}.cutadapt.txt"
        else:
            logger.debug(f"Output name override: {self.output_name_override}")
            self.output_file = self.output_dir / f"{self.output_name_override}.cutadapt.fastq"
            self.json_report = self.output_dir / f"{self.output_name_override}.cutadapt.json"
            self.txt_report = self.output_dir / f"{self.output_name_override}.cutadapt.txt"

        # This needs to be last so the output files are set correctly!
        if subsample > 0:
            max_reads = subsample
            tmp_file = NamedTemporaryFile(delete=False, suffix=".fastq")
            logger.info(f"Subsampling {max_reads} reads from input.")
            with open(self.read, "r") as r, open(tmp_file.name, "w") as out:
                for i in range(max_reads):
                    for _ in range(4):
                        out.write(r.readline())
            self.read = Path(tmp_file.name)
        self.double_tagged_reads = -1
        self.single_tagged_reads = -1
        self.untagged_reads = -1
        self.frac_tagged_reads = -1

    def run(self, threads=1, extra_params=None, quiet=False, **kwargs) -> Tuple[str, str]:
        cutadapt_call = [
            executable, "-m", "cutadapt",
            f"--action={self.action}",
            "-j", str(threads),
        ]
        # For single-end read
        cutadapt_call.extend(["-g", f"{self.adapter_name}={self.adapter_sequence}"])
        if self.adapter_sequence2 is not None:
            cutadapt_call.extend(["-g", f"{self.adapter_name2}={self.adapter_sequence2}"])
        add_to_cutadapt_call = [
            "--overlap", str(self.min_overlap),
            "--error-rate", str(self.error_rate),
            "--rename", "{id} Barcode={adapter_name} Seq={match_sequence} {comment}",
            "--minimum-length", "1",  # Discard reads with zero length after trimming
            str(self.read),
            "-o", str(self.output_file),
            f"--json={self.json_report}",
        ]
        cutadapt_call.extend(add_to_cutadapt_call)
        if self.times != 1:
            cutadapt_call.extend(["--times", str(self.times)])
        if self.reverse_complement_search:
            cutadapt_call.append("--rc")
        if isinstance(extra_params, list) or isinstance(extra_params, tuple):
            cutadapt_call.extend(extra_params)
        elif isinstance(extra_params, str):
            cutadapt_call.append(extra_params)
        if quiet:
            cutadapt_call.append("--quiet")
        logger.debug(f"Running cutadapt with the following command: {' '.join(cutadapt_call)}")
        result = run(cutadapt_call, capture_output=True, text=True)
        logger.debug(f"Cutadapt completed with return code {result.returncode}. ")
        if result.returncode != 0:
            raise RuntimeError(f"Cutadapt failed: {result.stderr}")
        with open(self.txt_report, "w") as f:
            f.write(result.stdout)
        return result.stdout, result.stderr

    def read_json_report(self) -> dict:
        if not self.json_report.exists():
            raise FileNotFoundError(f"JSON report {self.json_report} does not exist.")
        with open(self.json_report, 'r') as f:
            return load(f)

    def print_quick_report(self):
        if not self.json_report.exists():
            raise ValueError("JSON report not set. Run the cutadapt command first. Or use read_json_report() to load it.")
        report = self.read_json_report()
        read_counts = report.get('read_counts', {})
        total_out_reads = read_counts.get('output', 0)
        total_found_adapters = read_counts.get('with_adapters', 0)
        quick_report = (
            f"Quick Report:\n"
            f"\tInput file: {self.read}\n"
            f"\tOutput file: {self.output_file}\n"
            f"\tAdapter sequence: {self.adapter_sequence}\n"
            f"\tTotal output reads: {total_out_reads:>9} reads\n"
            f"\tTotal reads with adapter: {total_found_adapters:>9} "
            f"reads ({total_found_adapters / total_out_reads:.2%})\n"
            f"\tJSON report: {self.json_report}\n"
            f"\tText report: {self.txt_report}\n"
        )
        logger.info(quick_report)

    def provide_mini_results_dict(self):
        """
        Provides a mini results dictionary with basic information about the run.
        """
        if not self.json_report.exists():
            raise ValueError("JSON report not set. Run the cutadapt command first. Or use read_json_report() to load it.")
        report = self.read_json_report()
        read_counts = report.get('read_counts', {})
        total_out_reads = read_counts.get('output', 0)
        total_found_adapters = read_counts.get('with_adapters', 0)

        return {
            "output_file": str(self.output_file),
            "adapter_sequence": self.adapter_sequence,
            "total_output_reads": total_out_reads,
            "total_found_adapters": total_found_adapters,
            "percentage_found": (total_found_adapters / total_out_reads * 100) if total_out_reads > 0 else 0,
            "overall_percentage_barcoded": self.frac_tagged_reads,
            "number_tagged": self.single_tagged_reads + self.double_tagged_reads,
            "number_double_tagged": self.double_tagged_reads,
            "number_single_tagged": self.single_tagged_reads,
            "number_untagged": self.untagged_reads,
            "whole_output_json": report,
        }

    def update_json_report(self):
        # We want to load the existing json report, update it with the new information, and then write it back out.
        if not self.json_report.exists():
            raise ValueError("JSON report not set. Run the cutadapt command first. Or use read_json_report() to load it.")
        report = self.read_json_report()
        add_on_info = {
            "adapter_name": self.adapter_name,
            "adapter_sequence": self.adapter_sequence,
            "output_file": str(self.output_file),
            "number_double_tagged": self.double_tagged_reads,
            "number_single_tagged": self.single_tagged_reads,
            "number_untagged": self.untagged_reads,
            "overall_percentage_barcoded": self.frac_tagged_reads,
            "number_tagged": self.single_tagged_reads + self.double_tagged_reads,
        }
        report.update({"additional_info": add_on_info})
        with open(self.json_report, 'w') as f:
            import json
            json.dump(report, f, indent=2)
        logger.success(f"JSON report updated successfully.")

    def split_tagged_and_untagged_reads(self,
                                        name_for_hits: str = "tagged",
                                        name_for_misses: str = "untagged",
                                        ) -> None:
        """
        Splits the reads into tagged and untagged based on the presence of adapter sequence.
        """
        import gzip
        from mimetypes import guess_type
        from functools import partial

        encoding = guess_type(self.output_file)[1]  # uses file extension
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

        # Create new filename stems with barcoded info embedded in the main filename
        if self.output_name_override:
            base_name = self.output_name_override
            tagged_fastq = self.output_file.parent / f"{base_name}_{name_for_hits}.cutadapt.fastq"
            untagged_fastq = self.output_file.parent / f"{base_name}_{name_for_misses}.cutadapt.fastq"
        else:
            name_no_ext = self.read.name.replace(".fastq", "").replace(".fq", "")
            tagged_fastq = self.output_file.parent / f"{name_no_ext}_{name_for_hits}.cutadapt.fastq"
            untagged_fastq = self.output_file.parent / f"{name_no_ext}_{name_for_misses}.cutadapt.fastq"

        with _open(self.output_file) as handle:
            total_reads = sum(1 for _ in SeqIO.parse(handle, "fastq"))

        with open(tagged_fastq, 'w') as tagged, \
             open(untagged_fastq, 'w') as untagged, \
             _open(self.output_file) as handle:
            iterator = tqdm.tqdm(enumerate(SeqIO.parse(handle, "fastq")), total=total_reads)
            adapter_names = {self.adapter_name}
            if self.dual_adapters:
                adapter_names.add(self.adapter_name2)
            single_tags, double_tags, untags = 0, 0, 0
            for i, record in iterator:
                comments = record.description.split()
                # Extract barcode info from comments
                barcode_info = None
                for comment in comments:
                    if comment.startswith("Barcode="):
                        barcode_info = comment.split("=")[1]
                        break

                if barcode_info and barcode_info in adapter_names:
                    # This is a tagged read
                    SeqIO.write(record, tagged, "fastq")
                    single_tags += 1
                else:
                    # This is an untagged read
                    SeqIO.write(record, untagged, "fastq")
                    untags += 1

                if (i % 100 == 0 and i > 0) or i == total_reads - 1:
                    iterator.set_description(f"Splitting reads: "
                                             f"{name_for_hits.title()}: {single_tags} ({single_tags / (i+1):.1%}) "
                                             f"{name_for_misses.title()}: {untags} ({untags / (i+1):.1%})")

        logger.debug(f"Tagged reads written to {tagged_fastq}.")
        logger.debug(f"Untagged reads written to {untagged_fastq}.")
        logger.success(f"{name_for_hits.title()}: {single_tags} ({single_tags / total_reads:.1%}) "
                    f"{name_for_misses.title()}: {untags} ({untags / total_reads:.1%})")

        # Delete the original output file
        try:
            self.output_file.unlink(missing_ok=True)
        except Exception as e:
            logger.info(f"Could not delete original output file: {e}")
            logger.info("You may want to delete it manually if you don't need it anymore.")

        self.double_tagged_reads = double_tags
        self.single_tagged_reads = single_tags
        self.untagged_reads = untags
        self.frac_tagged_reads = (single_tags + double_tags) / (single_tags + double_tags + untags) if (single_tags + double_tags + untags) > 0 else 0


def print_version_info():
    """Print version information for the wrapper and its dependencies"""
    import subprocess
    import sys

    print("portable_cutadapt_wrapper.py version: 1.0.0")
    print(f"Python: {sys.version.split()[0]}")

    try:
        # Get cutadapt version
        result = subprocess.run([sys.executable, "-m", "cutadapt", "--version"],
                              capture_output=True, text=True)
        if result.returncode == 0:
            print(f"cutadapt: {result.stdout.strip()}")
        else:
            print("cutadapt: version unavailable")
    except Exception:
        print("cutadapt: version unavailable")

    # Get other dependency versions
    dependencies = ["pyyaml", "tqdm", "biopython", "loguru"]
    for dep in dependencies:
        try:
            import importlib.metadata
            version = importlib.metadata.version(dep)
            print(f"{dep}: {version}")
        except Exception:
            print(f"{dep}: version unavailable")

def _create_parser():
    """Create the argument parser (internal function)."""
    import argparse

    parser = argparse.ArgumentParser(description="Run cutadapt on single-end or paired-end FASTQ files.")
    parser.add_argument("fastq_r1", type=Path, help="Path to the FASTQ file (single-end) or first FASTQ file (R1, paired-end).")
    parser.add_argument("fastq_r2", type=Path, nargs='?', default=None, help="Path to the second FASTQ file (R2, paired-end only). If omitted, runs in single-end mode.")
    parser.add_argument("-c", "--config-file", type=Path, default=None,
                        help="Path to a YAML configuration file with default parameters.")
    parser.add_argument("--output-dir", type=Path, required=False, default=None,
                        help="Directory to save output files. Defaults to the directory of the input files.")
    parser.add_argument("--do-not-make-output-dir", action='store_true',
                        help="Do not create the output directory if it does not exist. "
                             "This will raise an error if the directory does not exist.")
    parser.add_argument("--output-name-override", type=str, default=None,
                        help="Override for output file names. [default, input_file_name+.cutadapt.<suffix>].")
    parser.add_argument("--adapter-name", type=str, default=None,
                        help="Name of the adapter sequence.")
    parser.add_argument("--adapter-sequence", type=str, required=False, default=None,
                        help="Adapter sequence to search for. Can be complete, or partial if you choose to use"
                             "the --anchor-sequence, --max-dist-from-read-start, and --convert options.")
    parser.add_argument("--adapter-name2", type=str, default=None,
                        help="Name of the second adapter sequence. Useful for dual-adapter scenarios.")
    parser.add_argument("--adapter-sequence2", type=str, default=None,
                        help="Second adapter sequence to search for. Useful for dual-adapter scenarios.")
    parser.add_argument("--anchor-sequence", type=str, default=None,
                        help="Anchor sequence to use for adapter matching. Default is None.")
    parser.add_argument("--max-dist-from-read-start", type=int,
                        default=None,
                        help="Maximum distance from the start of the read to search for the adapter sequence. "
                             "Default is None, which means.")
    parser.add_argument("--convert-CtoY", action='store_true',
                        help="Convert C to Y in the adapter sequence. Important for EM treated libraries.")
    parser.add_argument("--convert-GtoR", action='store_true',
                        help="Convert G to R in the adapter sequence. May be important for EM treated libraries?")
    parser.add_argument("--error-rate", type=float, default=0.0,
                        help="Error rate for adapter matching.")
    parser.add_argument("--reverse-complement-search", action='store_true',
                        help="Search for reverse complement of adapter.")
    parser.add_argument("--min-overlap", type=int, default=10,
                        help="Minimum overlap length for adapter matching.")
    parser.add_argument("--action", type=str,
                        choices=["trim", "retain", "none", "mask", "crop", "lowercase"],
                        default="trim",
                        help="Action to perform on reads with adapters.")
    parser.add_argument("--flip-reads", action='store_true',
                        help="Flip R1 and R2 reads.")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads to use.")
    parser.add_argument("--suffix-for-tagged-read-files", type=str, default="tagged",
                        help="Suffix for tagged read files. Default is 'tagged'.")
    parser.add_argument("--suffix-for-untagged-read-files", type=str, default="untagged",
                        help="Suffix for untagged read files. Default is 'untagged'.")
    parser.add_argument("--extra-params", nargs='*', default=None,
                        help="Additional parameters to pass to cutadapt.")
    parser.add_argument("--do-not-split", action='store_true',
                        help="Do not split reads into tagged and untagged. "
                             "This will keep the output as a single file.")
    parser.add_argument("--subsample", type=int, default=-1,
                        help="Subsample the input files to this many pairs. "
                             "Useful for testing with smaller datasets. Default is -1 (no subsampling).")
    parser.add_argument("--allow-indels", action='store_true',
                        help="Allow indels in the adapter sequence.")
    parser.add_argument("--log-path", type=Path, default=None,
                        help="Path to the logging file. Defaults to not writing log to a file.")
    parser.add_argument(
        "--log-level",
        "-l",
        type=str.upper,  # Convert input to uppercase for case-insensitivity
        choices=VALID_LOG_LEVELS,
        default="INFO",  # Default log level
        help=f"Set the logging level. Choices: {', '.join(VALID_LOG_LEVELS)} (default: INFO)",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Show version information for the wrapper and its dependencies"
    )
    return parser

def parse_args():
    import sys

    # Check for version flag first, before setting up full argument parsing
    if '--version' in sys.argv:
        print_version_info()
        sys.exit(0)

    parser = _create_parser()
    args = parser.parse_args()

    # Remove default handler and add a new one with the specified level
    logger.remove()
    logger.add(sys.stderr, level=args.log_level)
    if args.log_path is not None and args.log_path.parent.is_dir():
        logger.add(args.log_path, level="DEBUG")
    return args

def parse_yaml_config(config_file: Path):
    import yaml
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config

def merge_config_with_args(args, config):
    """
    Merges command line arguments with configuration file parameters.
    Precedence: Command line args > Config file params > Argparse defaults
    """
    # Get the parser and extract default values programmatically
    parser = _create_parser()
    defaults = parser.parse_args(['dummy1'])  # Parse with dummy positional arg (R2 is optional)
    default_values = vars(defaults)

    # from pprint import pprint
    merged_args = vars(args).copy()
    # pprint("Command line args:", merged_args)
    # pprint("Config file:", config)
    # pprint("Defaults:", default_values)

    for key, config_value in config.items():
        # Only override with config value if:
        # 1. The key exists in merged_args, AND
        # 2. The current value equals the default value (meaning user didn't specify it)
        if key in merged_args and merged_args[key] == default_values.get(key):
            merged_args[key] = config_value
            # logger.debug(f"Using config value for {key}: {config_value}")

    # pprint("Final merged args:", merged_args)
    return merged_args

def main(test_dict_override=None):
    if test_dict_override is not None:
        # For testing purposes, we can override the args with a dictionary
        merged_args = test_dict_override
    else:
        args = parse_args()
        if args.config_file:
            config = parse_yaml_config(args.config_file)
            merged_args = merge_config_with_args(args, config)
        else:
            merged_args = vars(args)

    # Check if running in single-end or paired-end mode
    is_single_end = merged_args["fastq_r2"] is None

    assert merged_args["fastq_r1"].exists(), f"FASTQ file {merged_args['fastq_r1']} does not exist."
    if not is_single_end:
        assert merged_args["fastq_r2"].exists(), f"R2 FASTQ file {merged_args['fastq_r2']} does not exist."
    assert merged_args["adapter_sequence"] is not None, "Adapter sequence must be provided in CLI or yaml."

    merged_args["adapter_sequence"] = assemble_adapter_sequence(
        merged_args["adapter_sequence"],
        anchor_sequence=merged_args.get("anchor_sequence", ""),
        n_length=merged_args.get("max_dist_from_read_start", 0),
        convert_CtoY=merged_args.get("convert_CtoY", False),
        convert_GtoR=merged_args.get("convert_GtoR", False)
    )
    if merged_args.get("adapter_sequence2", None) is not None:
        merged_args["adapter_sequence2"] = assemble_adapter_sequence(
            merged_args["adapter_sequence2"],
            anchor_sequence=merged_args.get("anchor_sequence", ""),
            n_length=merged_args.get("max_dist_from_read_start", 0),
            convert_CtoY=merged_args.get("convert_CtoY", False),
            convert_GtoR=merged_args.get("convert_GtoR", False)
        )

    # Initialize the appropriate runner based on mode
    if is_single_end:
        logger.info(f"Initializing SingleEndCutAdaptRunner (single-end mode)")
        logger.debug(merged_args)
        # Rename fastq_r1 to fastq for SingleEndCutAdaptRunner
        merged_args["fastq"] = merged_args.pop("fastq_r1")
        merged_args.pop("fastq_r2", None)  # Remove fastq_r2 if present
        merged_args.pop("flip_reads", None)  # Remove flip_reads (not applicable to single-end)
        runner = SingleEndCutAdaptRunner(**merged_args)
    else:
        logger.info(f"Initializing PairedEndCutAdaptRunner (paired-end mode)")
        logger.debug(merged_args)
        runner = PairedEndCutAdaptRunner(**merged_args)

    mode_name = "SingleEndCutAdaptRunner" if is_single_end else "PairedEndCutAdaptRunner"
    logger.info(f"Starting {mode_name}. Searching for the following adapter sequence(s):")
    if merged_args.get("adapter_sequence2", None) is not None:
        logger.info(f"1: {runner.adapter_sequence}")
        logger.info(f"2: {runner.adapter_sequence2}")
    else:
        logger.info(runner.adapter_sequence)
    result_stdout, result_stderr = runner.run(**merged_args)
    runner.print_quick_report()
    if not merged_args.get("do_not_split", False):
        logger.info("Splitting reads into tagged and untagged based on adapter presence.")
        runner.split_tagged_and_untagged_reads(name_for_hits=merged_args['suffix_for_tagged_read_files'],
                                               name_for_misses=merged_args['suffix_for_untagged_read_files'])
        runner.update_json_report()
    else:
        logger.info("Not splitting reads into tagged and untagged. Keeping the output as a single file.")

def assemble_adapter_sequence(adapter_base_seq: str,
                              anchor_sequence: str = "",
                              n_length: int = 0,
                              convert_CtoY: bool = False,
                              convert_GtoR: bool = False) -> str:
    """
    Assembles an adapter sequence based on the provided parameters.
    :param adapter_base_seq:
    :param anchor_sequence:
    :param n_length:
    :param CtoY:
    :param GtoR:
    :return:
    """
    if convert_CtoY:
        adapter_base_seq = adapter_base_seq.replace("C", "Y")
    if convert_GtoR:
        adapter_base_seq = adapter_base_seq.replace("G", "R")
    if n_length > 0:
        n_sequence = "N{" + str(n_length) + "}" if n_length > 1 else "N"
    else:
        n_sequence = ""
    if anchor_sequence == "" or not anchor_sequence:
        adapter_sequence = adapter_base_seq
    else:
        adapter_sequence = f"{anchor_sequence}{n_sequence}{adapter_base_seq}"
    return adapter_sequence

if __name__ == "__main__":
    # Test with:
    # ./portable_cutadapt_wrapper.py test_data/mot23/Ca3-RNA-NNSR_09Z_S33_L001_R1_001.first400k.fastq.gz test_data/mot23/Ca3-RNA-NNSR_09Z_S33_L001_R2_001.first400k.fastq.gz -c runner_settings_readStart_EM_e0.0.yaml --output-dir ./test_data/mot23/portable_cutadapt
    main()
