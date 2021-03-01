#!/bin/bash
#Parâmetros
# ./mdf 25 1000 500 250 253 0
#        |   |   |   |   |  |-> 1 grava em arquivo/0 não grava saída no arquivo
#        |   |   |   |   |----> temperatura máxima
#        |   |---|---|--------> tamanho do domínio 3d (exemplo: 1000 pontos em x, 500 em y e 250 em z)
#        |--------------------> passos de tempo
time ./mdf 25 1000 500 250 253 0
