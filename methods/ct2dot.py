# SPDX-License-Identifier: GPL-2.0-only
#
# Python port/adaptation of pseudoknot-aware dot-bracket conversion
# functionality from RNAstructure.
#
# Based on functionality implemented in RNAstructure's structure.cpp,
# including pseudoknot detection, ranking, non-crossing pair selection,
# and extended dot-bracket assignment.
#
# Python adaptation and extensions:
# Copyright (C) 2026 Jingwen Liu
#
# Modifications and extensions include:
# - BPSEQ input support
# - support for multiple CT header formats
# - batch conversion
# - logging and argparse CLI
# - support for long RNA structures
#
# Modified: 2026

"""
ct2dot.py - Convert CT (Connectivity Table) and BPSEQ files to dot-bracket notation

This Python implementation provides pseudoknot-aware dot-bracket conversion
based on functionality from RNAstructure, with additional support for BPSEQ
input, multiple CT header formats, batch conversion, and long RNA structures.

CT File Format:
- Line 1: <sequence_length> <structure_name>
- Subsequent lines: <index> <base> <prev_index> <next_index> <pair_index> <natural_index>
  - pair_index = 0 means unpaired
  - pair_index > 0 means paired with that position

BPSEQ File Format:
- Each line: <index> <base> <pair_index>
  - pair_index = 0 means unpaired
  - pair_index > 0 means paired with that position

Dot-Bracket Notation:
- '.' = unpaired base
- '(' = paired base (5' side, i < j where j is the pairing partner)
- ')' = paired base (3' side, i > j where j is the pairing partner)


"""

import os
import sys
import argparse
from pathlib import Path
import logging
import re


class BPSEQStructure:
    """Class to represent and manipulate BPSEQ structure data."""
    
    def __init__(self):
        self.length = 0
        self.name = ""
        self.sequence = []
        self.pairs = []  # pairs[i] = j means position i pairs with position j (0 = unpaired)
    
    def read_bpseq_file(self, bpseq_file_path):
        """
        Read a BPSEQ file and parse its contents.
        
        BPSEQ Format:
        - Each line: <index> <base> <pair_index>
        - pair_index = 0 means unpaired
        - pair_index > 0 means paired with that position
        
        Args:
            bpseq_file_path: Path to the BPSEQ file
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            with open(bpseq_file_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) == 0:
                logging.error(f"Empty BPSEQ file: {bpseq_file_path}")
                return False
            
            # Parse BPSEQ format
            entries = []
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) < 3:
                    logging.warning(f"Skipping invalid line: {line}")
                    continue
                
                try:
                    index = int(parts[0])
                    base = parts[1]
                    pair_index = int(parts[2])
                    entries.append((index, base, pair_index))
                except ValueError as e:
                    logging.warning(f"Skipping line with parse error: {line} - {e}")
                    continue
            
            if not entries:
                logging.error(f"No valid entries found in BPSEQ file: {bpseq_file_path}")
                return False
            
            self.length = len(entries)
            self.name = Path(bpseq_file_path).stem
            
            # Initialize arrays (1-indexed)
            self.sequence = [''] * (self.length + 1)
            self.pairs = [0] * (self.length + 1)
            
            # Fill arrays
            for index, base, pair_index in entries:
                if 1 <= index <= self.length:
                    self.sequence[index] = base
                    self.pairs[index] = pair_index
                else:
                    logging.warning(f"Index {index} out of range [1, {self.length}]")
            
            logging.info(f"Successfully read BPSEQ file: {bpseq_file_path}")
            logging.info(f"  Name: {self.name}")
            logging.info(f"  Length: {self.length} nt")
            
            return True
            
        except Exception as e:
            logging.error(f"Error reading BPSEQ file {bpseq_file_path}: {str(e)}")
            return False
    
    def to_dot_bracket(self, use_extended_notation=True):
        """
        Convert the BPSEQ structure to dot-bracket notation.
        
        Args:
            use_extended_notation: If True, uses extended notation for pseudoknots
                                   with multiple bracket types: ()<>{}[]AaBbCc...
                                   If False, uses only () for all pairs
        
        Returns:
            str: Dot-bracket notation string
        """
        if self.length == 0:
            return ""
        
        if not use_extended_notation:
            return self._to_simple_dot_bracket()
        
        return self._to_extended_dot_bracket_rnastructure()
    
    def _to_simple_dot_bracket(self):
        """Convert to simple dot-bracket using only () for all pairs."""
        dot_bracket = ['.'] * (self.length + 1)
        
        for i in range(1, self.length + 1):
            pair_j = self.pairs[i]
            
            if pair_j == 0:
                dot_bracket[i] = '.'
            elif i < pair_j:
                dot_bracket[i] = '('
            elif i > pair_j:
                dot_bracket[i] = ')'
            else:
                logging.warning(f"Invalid self-pairing at position {i}")
                dot_bracket[i] = '.'
        
        return ''.join(dot_bracket[1:])
    
    def _to_extended_dot_bracket_rnastructure(self):
        """
        Convert to extended dot-bracket notation using RNAstructure's algorithm.
        Bracket symbols: ()<>{}[]AaBbCcDd...
        """
        brackets = "()<>{}[]AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz"
        max_level = len(brackets) // 2
        
        ranks = self._get_pseudoknot_ranks()
        
        dot_bracket = ['.'] * (self.length + 1)
        
        for i in range(1, self.length + 1):
            pair_j = self.pairs[i]
            level = min(ranks[i], max_level) - 1
            
            if pair_j > i:
                dot_bracket[i] = brackets[2 * level]
            elif pair_j == 0:
                dot_bracket[i] = '.'
            else:
                dot_bracket[i] = brackets[2 * level + 1]
        
        return ''.join(dot_bracket[1:])
    
    def _get_pseudoknot_ranks(self):
        """Calculate pseudoknot rank for each position."""
        ranks = [0 if self.pairs[i] == 0 else 1 for i in range(len(self.pairs))]
        
        current_pairs = self.pairs[:]
        
        while self._has_pseudoknots(current_pairs):
            pseudoknots = self._find_pseudoknots(current_pairs)
            
            for i in range(len(ranks)):
                if pseudoknots[i] != 0:
                    ranks[i] += 1
            
            current_pairs = pseudoknots
        
        return ranks
    
    def _has_pseudoknots(self, pairs):
        """Check if the given pairing has any pseudoknots (crossing bonds)."""
        length = len(pairs)
        stack = [(1, length - 1)]
        
        while stack:
            i, j = stack.pop()
            
            while i <= j and pairs[i] == 0:
                i += 1
            
            if i > j:
                continue
            
            k = pairs[i]
            if k < i:
                logging.error(f"Logic error: 5' end encountered at {i}")
                continue
            
            if k > j:
                return True
            
            if k + 1 <= j:
                stack.append((k + 1, j))
            if i + 1 <= k - 1:
                stack.append((i + 1, k - 1))
        
        return False
    
    def _find_pseudoknots(self, pairs):
        """Find pseudoknots using dynamic programming."""
        length = len(pairs)
        if length == 0:
            return [0] * length
        
        outer = {}
        trace = {}
        
        for i in range(1, length):
            outer[(i, i)] = 0
        
        for n in range(1, length):
            for i in range(1, length - n):
                j = i + n
                outer[(i, j)] = outer.get((i + 1, j), 0)
                trace[(i, j)] = False
                
                k = pairs[i]
                if k != 0 and k > i and k <= j:
                    tmp = 1
                    if i + 1 <= k - 1:
                        tmp += outer.get((i + 1, k - 1), 0)
                    if k + 1 <= j:
                        tmp += outer.get((k + 1, j), 0)
                    
                    if tmp >= outer[(i, j)]:
                        outer[(i, j)] = tmp
                        trace[(i, j)] = True
        
        results = pairs[:]
        stack = [(1, length - 1)]
        
        while stack:
            i, j = stack.pop()
            
            while i < j and not trace.get((i, j), False):
                i += 1
            
            if i >= j:
                continue
            
            k = pairs[i]
            results[i] = -k
            results[k] = -pairs[k]
            
            if i + 1 < k - 1:
                stack.append((i + 1, k - 1))
            if k + 1 < j:
                stack.append((k + 1, j))
        
        knotted = [0] * length
        for i in range(1, length):
            if results[i] > 0:
                knotted[i] = results[i]
        
        return knotted
    
    def get_sequence_string(self):
        """Get the sequence as a string."""
        return ''.join(self.sequence[1:])
    
    def validate_structure(self):
        """Validate the structure for consistency."""
        errors = []
        
        for i in range(1, self.length + 1):
            j = self.pairs[i]
            if j > 0:
                if j > self.length:
                    errors.append(f"Position {i} pairs with out-of-bounds position {j}")
                elif self.pairs[j] != i:
                    errors.append(f"Inconsistent pairing: {i}<->{j}, but {j}<->{self.pairs[j]}")
        
        is_valid = len(errors) == 0
        return is_valid, errors


class CTStructure:
    """Class to represent and manipulate CT (Connectivity Table) structure data."""
    
    def __init__(self):
        self.length = 0
        self.name = ""
        self.sequence = []
        self.pairs = []  # pairs[i] = j means position i pairs with position j (0 = unpaired)
    
    def read_ct_file(self, ct_file_path):
        """
        Read a CT file and parse its contents.
        Supports two header formats:
        1. Standard: <sequence_length> <structure_name>
        2. FASTA-like: >seq length: 890 seq name: IAV-vRNA8
        
        Args:
            ct_file_path: Path to the CT file
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            with open(ct_file_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) == 0:
                logging.error(f"Empty CT file: {ct_file_path}")
                return False
            
            # Parse header line - handle multiple formats
            header_line = lines[0].strip()
            line_index = 0
            
            # Check if header starts with '>' (FASTA-like format)
            if header_line.startswith('>'):
                # Format: >seq length: 890 seq name: IAV-vRNA8
                self.name = header_line[1:]  # Remove '>'
                line_index = 1
                
                # Try to extract length from the header
                # Look for "length: <number>" pattern
                match = re.search(r'length:\s*(\d+)', header_line)
                if match:
                    self.length = int(match.group(1))
                else:
                    # If not found, will infer from data lines
                    self.length = 0
            else:
                # Standard format: <sequence_length> <structure_name>
                header = header_line.split()
                if len(header) < 1:
                    logging.error(f"Invalid header in CT file: {ct_file_path}")
                    return False
                
                try:
                    self.length = int(header[0])
                    self.name = " ".join(header[1:]) if len(header) > 1 else "Unknown"
                except ValueError:
                    logging.error(f"Invalid header format in CT file: {ct_file_path}")
                    return False
                
                line_index = 1
            
            # Initialize arrays with provisional length
            # We'll resize if needed
            if self.length > 0:
                self.sequence = [''] * (self.length + 1)  # 1-indexed
                self.pairs = [0] * (self.length + 1)       # 1-indexed
            else:
                self.sequence = ['']
                self.pairs = [0]
            
            # Parse structure lines
            max_index = 0
            for i in range(line_index, len(lines)):
                line = lines[i].strip()
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) < 5:
                    logging.warning(f"Skipping invalid line {i}: {line}")
                    continue
                
                try:
                    index = int(parts[0])
                    base = parts[1]
                    pair_index = int(parts[4])
                    
                    max_index = max(max_index, index)
                    
                    # Resize arrays if necessary
                    if index >= len(self.sequence):
                        new_size = index + 1
                        self.sequence.extend([''] * (new_size - len(self.sequence)))
                        self.pairs.extend([0] * (new_size - len(self.pairs)))
                    
                    self.sequence[index] = base
                    self.pairs[index] = pair_index
                except (ValueError, IndexError) as e:
                    logging.warning(f"Error parsing line {i}: {line} - {e}")
                    continue
            
            # Update length if it wasn't in header
            if self.length == 0:
                self.length = max_index
            
            logging.info(f"Successfully read CT file: {ct_file_path}")
            logging.info(f"  Structure: {self.name}")
            logging.info(f"  Length: {self.length} nt")
            
            return True
            
        except Exception as e:
            logging.error(f"Error reading CT file {ct_file_path}: {str(e)}")
            return False
    
    def to_dot_bracket(self, use_extended_notation=True):
        """
        Convert the CT structure to dot-bracket notation.
        Following RNAstructure's algorithm for pseudoknot handling.
        
        Args:
            use_extended_notation: If True, uses extended notation for pseudoknots
                                   with multiple bracket types: ()<>{}[]AaBbCc...
                                   If False, uses only () for all pairs
        
        Returns:
            str: Dot-bracket notation string
        """
        if self.length == 0:
            return ""
        
        if not use_extended_notation:
            # Simple mode: use only () for all pairs
            return self._to_simple_dot_bracket()
        
        # Extended mode: detect pseudoknots and use multiple bracket types
        # This follows RNAstructure's algorithm
        return self._to_extended_dot_bracket_rnastructure()
    
    def _to_simple_dot_bracket(self):
        """Convert to simple dot-bracket using only () for all pairs."""
        dot_bracket = ['.'] * (self.length + 1)
        
        for i in range(1, self.length + 1):
            pair_j = self.pairs[i]
            
            if pair_j == 0:
                dot_bracket[i] = '.'
            elif i < pair_j:
                dot_bracket[i] = '('
            elif i > pair_j:
                dot_bracket[i] = ')'
            else:
                logging.warning(f"Invalid self-pairing at position {i}")
                dot_bracket[i] = '.'
        
        return ''.join(dot_bracket[1:])
    
    def _to_extended_dot_bracket_rnastructure(self):
        """
        Convert to extended dot-bracket notation using RNAstructure's algorithm.
        Bracket symbols: ()<>{}[]AaBbCcDd...
        """
        # RNAstructure's bracket symbols
        brackets = "()<>{}[]AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz"
        max_level = len(brackets) // 2
        
        # Get pseudoknot ranks for each position
        ranks = self._get_pseudoknot_ranks()
        
        # Build dot-bracket string
        dot_bracket = ['.'] * (self.length + 1)
        
        for i in range(1, self.length + 1):
            pair_j = self.pairs[i]
            # Use min to cap at max_level
            level = min(ranks[i], max_level) - 1
            
            if pair_j > i:
                # Opening bracket
                dot_bracket[i] = brackets[2 * level]
            elif pair_j == 0:
                # Unpaired
                dot_bracket[i] = '.'
            else:
                # Closing bracket (pair_j < i)
                dot_bracket[i] = brackets[2 * level + 1]
        
        return ''.join(dot_bracket[1:])
    
    def _get_pseudoknot_ranks(self):
        """
        Calculate pseudoknot rank for each position.
        Rank 0 = unpaired
        Rank 1 = non-crossing pair
        Rank 2 = first-order pseudoknot
        Rank 3 = second-order pseudoknot, etc.
        
        This implements RNAstructure's GetPseudoknotRanks algorithm.
        """
        # Initialize ranks: 0 for unpaired, 1 for all pairs
        ranks = [0 if self.pairs[i] == 0 else 1 for i in range(len(self.pairs))]
        
        # Work with a copy of pairs
        current_pairs = self.pairs[:]
        
        # Iteratively find pseudoknots and increment ranks
        while self._has_pseudoknots(current_pairs):
            # Find pseudoknots in current_pairs
            pseudoknots = self._find_pseudoknots(current_pairs)
            
            # Increment rank for each pseudoknot position
            for i in range(len(ranks)):
                if pseudoknots[i] != 0:
                    ranks[i] += 1
            
            # Update current_pairs to only contain pseudoknots for next iteration
            current_pairs = pseudoknots
        
        return ranks
    
    def _has_pseudoknots(self, pairs):
        """
        Check if the given pairing has any pseudoknots (crossing bonds).
        Implements RNAstructure's hasPseudoknots algorithm.
        """
        length = len(pairs)
        stack = [(1, length - 1)]
        
        while stack:
            i, j = stack.pop()
            
            # Find next paired position
            while i <= j and pairs[i] == 0:
                i += 1
            
            if i > j:
                continue
            
            k = pairs[i]
            if k < i:
                logging.error(f"Logic error: 5' end encountered at {i}")
                continue
            
            # If 3' end is outside interval, it's a crossing bond
            if k > j:
                return True
            
            # Push intervals for further processing
            if k + 1 <= j:
                stack.append((k + 1, j))
            if i + 1 <= k - 1:
                stack.append((i + 1, k - 1))
        
        return False
    
    def _find_pseudoknots(self, pairs):
        """
        Find pseudoknots using dynamic programming.
        Returns two sets: optimal (non-crossing) and knotted (crossing).
        We only need the knotted set here.
        
        Implements RNAstructure's findPseudoknots algorithm.
        """
        length = len(pairs)
        if length == 0:
            return [0] * length
        
        # Dynamic programming table
        # outer[i][j] = maximum number of non-crossing bonds in interval [i,j]
        outer = {}
        trace = {}
        
        # Initialize diagonal
        for i in range(1, length):
            outer[(i, i)] = 0
        
        # Fill DP table
        for n in range(1, length):
            for i in range(1, length - n):
                j = i + n
                outer[(i, j)] = outer.get((i + 1, j), 0)
                trace[(i, j)] = False
                
                k = pairs[i]
                if k != 0 and k > i and k <= j:
                    tmp = 1
                    if i + 1 <= k - 1:
                        tmp += outer.get((i + 1, k - 1), 0)
                    if k + 1 <= j:
                        tmp += outer.get((k + 1, j), 0)
                    
                    if tmp >= outer[(i, j)]:
                        outer[(i, j)] = tmp
                        trace[(i, j)] = True
        
        # Backtrace to find optimal set
        results = pairs[:]
        stack = [(1, length - 1)]
        
        while stack:
            i, j = stack.pop()
            
            while i < j and not trace.get((i, j), False):
                i += 1
            
            if i >= j:
                continue
            
            k = pairs[i]
            results[i] = -k  # Mark as in optimal set
            results[k] = -pairs[k]
            
            if i + 1 < k - 1:
                stack.append((i + 1, k - 1))
            if k + 1 < j:
                stack.append((k + 1, j))
        
        # Return knotted pairs (those NOT in optimal set)
        knotted = [0] * length
        for i in range(1, length):
            if results[i] > 0:  # Not in optimal set
                knotted[i] = results[i]
        
        return knotted
    
    def get_sequence_string(self):
        """
        Get the sequence as a string.
        
        Returns:
            str: Sequence string
        """
        return ''.join(self.sequence[1:])  # Skip index 0
    
    def validate_structure(self):
        """
        Validate the structure for consistency.
        
        Returns:
            tuple: (is_valid, error_messages)
        """
        errors = []
        
        # Check pairing consistency: if i pairs with j, then j should pair with i
        for i in range(1, self.length + 1):
            j = self.pairs[i]
            if j > 0:
                if j > self.length:
                    errors.append(f"Position {i} pairs with out-of-bounds position {j}")
                elif self.pairs[j] != i:
                    errors.append(f"Inconsistent pairing: {i}<->{j}, but {j}<->{self.pairs[j]}")
        
        is_valid = len(errors) == 0
        return is_valid, errors


def convert_ct_to_dot(ct_file_path, output_file_path, format_type='full', structure_number=1, use_extended_notation=True):
    """
    Convert a CT file to dot-bracket format.
    
    Args:
        ct_file_path: Path to input CT file
        output_file_path: Path to output dot-bracket file
        format_type: Output format ('simple', 'side', 'multi', 'full')
        structure_number: Which structure to convert (1-indexed, -1 for all)
        use_extended_notation: Use extended notation for pseudoknots (()[]{}<>)
        
    Returns:
        int: 0 on success, 1 on failure
    """
    # For now, we only support single structure CT files (structure_number parameter ignored)
    # The CT files in hasps directory contain single structures
    
    structure = CTStructure()
    
    if not structure.read_ct_file(ct_file_path):
        return 1
    
    # Validate structure
    is_valid, errors = structure.validate_structure()
    if not is_valid:
        logging.warning(f"Structure validation warnings for {ct_file_path}:")
        for error in errors:
            logging.warning(f"  {error}")
    
    # Convert to dot-bracket
    dot_bracket = structure.to_dot_bracket(use_extended_notation=use_extended_notation)
    sequence = structure.get_sequence_string()
    
    # Write output based on format type
    try:
        with open(output_file_path, 'w') as f:
            if format_type == 'simple':
                # Only structure line
                f.write(f"{dot_bracket}\n")
            
            elif format_type == 'side':
                # Structure with label on the side
                f.write(f"{dot_bracket}  {structure.name}\n")
            
            elif format_type == 'multi':
                # Title line + structure line
                f.write(f">{structure.name}\n")
                f.write(f"{dot_bracket}\n")
            
            else:  # 'full' is default
                # Title line + sequence line + structure line
                f.write(f">{structure.name}\n")
                f.write(f"{sequence}\n")
                f.write(f"{dot_bracket}\n")
        
        logging.info(f"Successfully wrote output to {output_file_path}")
        return 0
        
    except Exception as e:
        logging.error(f"Error writing output file {output_file_path}: {str(e)}")
        return 1


def convert_bpseq_to_dot(bpseq_file_path, output_file_path, format_type='full', use_extended_notation=True):
    """
    Convert a BPSEQ file to dot-bracket format.
    
    Args:
        bpseq_file_path: Path to input BPSEQ file
        output_file_path: Path to output dot-bracket file
        format_type: Output format ('simple', 'side', 'multi', 'full')
        use_extended_notation: Use extended notation for pseudoknots (()[]{}<>)
        
    Returns:
        int: 0 on success, 1 on failure
    """
    structure = BPSEQStructure()
    
    if not structure.read_bpseq_file(bpseq_file_path):
        return 1
    
    # Validate structure
    is_valid, errors = structure.validate_structure()
    if not is_valid:
        logging.warning(f"Structure validation warnings for {bpseq_file_path}:")
        for error in errors:
            logging.warning(f"  {error}")
    
    # Convert to dot-bracket
    dot_bracket = structure.to_dot_bracket(use_extended_notation=use_extended_notation)
    sequence = structure.get_sequence_string()
    
    # Write output based on format type
    try:
        with open(output_file_path, 'w') as f:
            if format_type == 'simple':
                # Only structure line
                f.write(f"{dot_bracket}\n")
            
            elif format_type == 'side':
                # Structure with label on the side
                f.write(f"{dot_bracket}  {structure.name}\n")
            
            elif format_type == 'multi':
                # Title line + structure line
                f.write(f">{structure.name}\n")
                f.write(f"{dot_bracket}\n")
            
            else:  # 'full' is default
                # Title line + sequence line + structure line
                f.write(f">{structure.name}\n")
                f.write(f"{sequence}\n")
                f.write(f"{dot_bracket}\n")
        
        logging.info(f"Successfully wrote output to {output_file_path}")
        return 0
        
    except Exception as e:
        logging.error(f"Error writing output file {output_file_path}: {str(e)}")
        return 1


def main():
    """Main function to handle command line interface."""
    
    parser = argparse.ArgumentParser(
        description='Convert CT/BPSEQ files to dot-bracket notation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input Format Types:
  ct      - CT (Connectivity Table) format (default)
  bpseq   - BPSEQ format

Output Format Types:
  simple  - Only structure line (no title or sequence)
  side    - Structure with label appended on the right
  multi   - Title line + structure line
  full    - Title line + sequence line + structure line (default)

Examples:
  %(prog)s input.ct output.dot
  %(prog)s input.ct output.dot -f simple
  %(prog)s input.bpseq output.dot --input-type bpseq
  %(prog)s input.bpseq output.dot --input-type bpseq -f full
  %(prog)s --batch-convert-dir hasps --output-dir psdot
  %(prog)s --batch-convert-dir bpseq_files --output-dir dot_output --input-type bpseq
        """
    )
    
    parser.add_argument('ct_file', nargs='?', help='Input CT/BPSEQ file path')
    parser.add_argument('output_file', nargs='?', help='Output dot-bracket file path')
    
    parser.add_argument('-f', '--format', 
                        choices=['simple', 'side', 'multi', 'full'],
                        default='full',
                        help='Output format type (default: full)')
    
    parser.add_argument('-t', '--input-type',
                        choices=['ct', 'bpseq'],
                        default='ct',
                        help='Input file format type (default: ct)')
    
    parser.add_argument('-n', '--structure-number',
                        type=int,
                        default=1,
                        help='Structure number to convert for CT files (1-indexed, -1 for all)')
    
    parser.add_argument('-b', '--batch-convert-dir',
                        help='Batch convert all files in this directory')
    
    parser.add_argument('-o', '--output-dir',
                        help='Output directory for batch conversion')
    
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Enable verbose logging')
    
    parser.add_argument('-q', '--quiet',
                        action='store_true',
                        help='Suppress unnecessary output')
    
    parser.add_argument('--simple-brackets',
                        action='store_true',
                        help='Use only () brackets (ignore pseudoknots)')
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.WARNING
    if args.verbose:
        log_level = logging.DEBUG
    elif not args.quiet:
        log_level = logging.INFO
    
    logging.basicConfig(
        level=log_level,
        format='%(levelname)s: %(message)s'
    )
    
    # Batch conversion mode
    if args.batch_convert_dir:
        input_dir = Path(args.batch_convert_dir)
        output_dir = Path(args.output_dir) if args.output_dir else input_dir / 'dot_output'
        
        if not input_dir.exists():
            logging.error(f"Input directory does not exist: {input_dir}")
            return 1
        
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Select file extension based on input type
        file_pattern = '*.bpseq' if args.input_type == 'bpseq' else '*.ct'
        input_files = list(input_dir.glob(file_pattern))
        
        if not input_files:
            logging.warning(f"No {file_pattern} files found in {input_dir}")
            return 0
        
        if not args.quiet:
            print(f"Converting {len(input_files)} {args.input_type.upper()} files from {input_dir} to {output_dir}")
        
        success_count = 0
        failed_count = 0
        
        for input_file in input_files:
            output_file = output_dir / f"{input_file.stem}.dot"
            
            if not args.quiet:
                print(f"  {input_file.name} -> {output_file.name}", end=' ... ')
            
            if args.input_type == 'bpseq':
                result = convert_bpseq_to_dot(
                    str(input_file),
                    str(output_file),
                    args.format,
                    use_extended_notation=not args.simple_brackets
                )
            else:
                result = convert_ct_to_dot(
                    str(input_file),
                    str(output_file),
                    args.format,
                    args.structure_number,
                    use_extended_notation=not args.simple_brackets
                )
            
            if result == 0:
                success_count += 1
                if not args.quiet:
                    print("✓")
            else:
                failed_count += 1
                if not args.quiet:
                    print("✗")
        
        if not args.quiet:
            print(f"\nConversion complete: {success_count} succeeded, {failed_count} failed")
        
        return 0 if failed_count == 0 else 1
    
    # Single file conversion mode
    else:
        if not args.ct_file or not args.output_file:
            parser.error("ct_file and output_file are required for single file conversion")
            return 1
        
        if not args.quiet:
            print(f"Converting {args.ct_file} to {args.output_file}")
        
        if args.input_type == 'bpseq':
            result = convert_bpseq_to_dot(
                args.ct_file,
                args.output_file,
                args.format,
                use_extended_notation=not args.simple_brackets
            )
        else:
            result = convert_ct_to_dot(
                args.ct_file,
                args.output_file,
                args.format,
                args.structure_number,
                use_extended_notation=not args.simple_brackets
            )
        
        if result == 0 and not args.quiet:
            print("Conversion complete.")
        
        return result


if __name__ == '__main__':
    sys.exit(main())
