function [y,fs] = AudioLoad_func(FILE)
params = load(FILE);
y = params.Aud.data;
fs = 48000;