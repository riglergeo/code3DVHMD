#!/bin/bash
modelo="$1"

extin=".in"
arq_in="$modelo$extin"

if [ -f "$arq_in" ]
then
  ./gera_poly $arq_in
else
  echo "Arquivo $arq_in não existe!"
  exit
fi

extpoly=".poly"
arq_poly="$modelo$extpoly"
tetgen -pq1.25 -keaA $arq_poly

./gera_arestas $modelo


