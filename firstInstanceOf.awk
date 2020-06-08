#!/bin/bash

awk 'BEGIN{id=""} (id!=$1) {print $0;id=$1}'
