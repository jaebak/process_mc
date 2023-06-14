#!/usr/bin/env python3
import subprocess
import os
import re
import string
import argparse
from multiprocessing import Pool
from datetime import datetime

def sluggify(text):
    # Remove leading/trailing whitespaces and convert to lowercase
    text = text.strip().lower()

    # Replace spaces and non-alphanumeric characters with hyphens
    text = re.sub(r"[^\w\s-]", "", text)
    text = re.sub(r"\s+", "-", text)

    # Remove any remaining non-ASCII characters
    text = "".join(char for char in text if char in string.ascii_letters or char.isdigit() or char == "-")

    return text

def run_command(command, print_output):
    # Create the log file name based on the sluggified command and current timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    sluggified_command = sluggify(command)
    log_file = os.path.join(args.log_dir, f"{sluggified_command}_{timestamp}.log")

    # Write the command at the top of the log file
    with open(log_file, "w", buffering=1) as file:
        file.write(f"Command: {command}\n")

        print(f'Running command : {command}\n  Log file: {log_file}')
        # Execute the command and capture the output in real-time
        with subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, bufsize=1) as process:
            # Write the output to the log file in real-time
            for line in process.stdout:
                file.write(line)
                if print_output:
                    print(line, end='', flush=True)
        print(f'Finished command : {command}\n  Log file: {log_file}')

def process_command(command):
    run_command(command, args.print_output)

if __name__ == "__main__":
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Run commands from a file and save the output to log files.")
    parser.add_argument("command_file", help="path to the file containing commands to run")
    parser.add_argument("--log-dir", default="command_logs", help="directory to save the log files")
    parser.add_argument("--print-output", action="store_true", help="print the output of the commands to the console")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Validate the log directory
    os.makedirs(args.log_dir, exist_ok=True)

    with open(args.command_file, "r") as file:
        commands = file.read().splitlines()

    # Remove lines starting with '#' to ignore them as comments
    commands = [command.strip() for command in commands if not command.startswith('#')]

    with Pool() as pool:
        pool.map(process_command, commands)
