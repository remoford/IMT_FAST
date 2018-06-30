#!/bin/bash
(pgrep icc_fast ; pgrep gcc_fast ) | xargs -n1 kill

