#!/bin/bash
#Syntax: commentOut <stringTo Search, e.g., version30x> <file>
#This will do a GLOBAL commenting out!

sed -i "/$1/ s:^://:g" $2
