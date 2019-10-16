#!/usr/bin/env bash

cd diff_cores_same_socket || exit
qsub osu_bw.job
qsub osu_bw.job
qsub osu_bw.job
qsub osu_latency.job
qsub osu_latency.job
qsub osu_latency.job
cd ../diff_nodes || exit
qsub osu_bw.job
qsub osu_bw.job
qsub osu_bw.job
qsub osu_latency.job
qsub osu_latency.job
qsub osu_latency.job
cd ../diff_socket_same_node || exit
qsub osu_bw.job
qsub osu_bw.job
qsub osu_bw.job
qsub osu_latency.job
qsub osu_latency.job
qsub osu_latency.job