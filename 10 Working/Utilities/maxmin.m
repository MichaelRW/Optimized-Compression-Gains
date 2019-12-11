function [ mm ] = maxmin( input )

mm = [max(max(max(input))) min(min(min(input)))];