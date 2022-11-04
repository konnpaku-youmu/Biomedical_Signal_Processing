clc;
clear;
close all;

load eegdata_artifacts.mat

eegplot_simple(eegdata, fs);
