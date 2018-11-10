#!/bin/bash

for outputFile in output/latest/*
do
	tail -f $outputFile &
done

