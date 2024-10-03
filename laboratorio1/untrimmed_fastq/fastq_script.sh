#!/bin/bash

for filename in *.fastq
	do ls $filename
	wc -l $filename
	echo Terminado
done

