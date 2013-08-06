#!/bin/bash

circos -conf etc/indels.conf -debug_group summary,timer,textplace > run.out
circos -conf etc/mobs.conf -debug_group summary,timer,textplace > run.out
circos -conf etc/mutations.conf -debug_group summary,timer,textplace > run.out
circos -conf etc/combined_circos.conf -debug_group summary,timer,textplace > run.out
