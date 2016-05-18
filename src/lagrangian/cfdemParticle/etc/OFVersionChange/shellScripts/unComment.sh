#!/bin/bash
#Syntax: unComment <stringTo Search, e.g., version30x> <file>
#This will do a GLOBAL uncomment!

sed -i "/$1/ s:^//::g" $2
