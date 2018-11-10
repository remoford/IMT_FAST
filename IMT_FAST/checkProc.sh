#!/bin/bash
#$ -t 1-40

hostname
uptime
grep "model name" /proc/cpuinfo | uniq -c
