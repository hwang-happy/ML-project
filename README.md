# ML-project
Proteins' Secondary Structure Prediction

Celem projektu jest przewidywanie struktury drugorzędowej białek na podstawie podanej sekwencji aminokwasowej. Struktury, które są brane pod uwagę to: helisy, kartki i skręty przypisane poszczególnym aminokwasom. \newline
Postanowiliśmy skupić się na implementacji kilku metod w celu ich porównania: algorytm Chou-Fasmana, ukryte Modele Markowa (HMM) oraz dwukierunkową sieć neuronową typu GRU. Wyniki metod są porównywane ze strukturami znajdującymi się w bazie DSSP.

Dokumentacja zawierająca opis użytych metod i parametrów znajduje się [tutaj](https://www.overleaf.com/read/wpqzxngmfgrn).

W celu łatwiejszego użycia programu został zaimplementowany również interfejs graficzny. Żeby go uruchomić wystarczy wpisać w linię komend polecenie 'pyhon GUI.py'.
