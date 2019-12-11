#!/bin/bash
#SBATCH -t 05-0:00
#SBATCH --mem=2000 
#SBATCH -A def-ibruce
octave test_speech_gains2a_65_ft.m
