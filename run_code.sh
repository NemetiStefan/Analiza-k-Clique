#!/bin/bash

# Compile the C program
gcc kClique.c -o program -lm
if [ $? -ne 0 ]; then
  echo "Compilation failed. Exiting."
  exit 1
fi

# Process 15 input files
for i in {1..15}; do
  input_file="input${i}.txt"
  output_file="output${i}.txt"

  # Check if the input file exists
  if [ ! -f "$input_file" ]; then
    echo "Input file $input_file not found. Skipping."
    continue
  fi

  # Run the program with the input file and save the output
  ./program < "$input_file" > "$output_file"

  if [ $? -eq 0 ]; then
    echo "Processed $input_file -> $output_file"
  else
    echo "Error processing $input_file"
  fi
done
