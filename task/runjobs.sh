#!/bin/bash

for file in jobs/*.job; do
    echo "$file"
    llsubmit "$file"
done
