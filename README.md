# Projekt z Nauczania Maszynowego - Proteins' Secondary Structure Prediction

Celem projektu jest przewidywanie struktury drugorzędowej białek na podstawie podanej sekwencji aminokwasowej. Struktury, które są brane pod uwagę to: helisy, kartki i skręty przypisane poszczególnym aminokwasom. \newline
Postanowiliśmy skupić się na implementacji kilku metod w celu ich porównania: algorytm Chou-Fasmana, ukryte Modele Markowa (HMM) oraz dwukierunkową sieć neuronową typu GRU. Wyniki metod są porównywane ze strukturami znajdującymi się w bazie DSSP.

Dokumentacja zawierająca opis użytych metod i parametrów znajduje się [tutaj](https://www.overleaf.com/read/wpqzxngmfgrn).

W celu łatwiejszego użycia programu został zaimplementowany również interfejs graficzny. Żeby go uruchomić wystarczy wpisać w linię komend polecenie `pyhon GUI.py`.

## Obsługa GUI
Dla poprawnego działania programu na wejście należy wkleić sekwencję aminokwasową w formacie FASTA (łącznie z nagłówkiem) lub wczytać ją z pliku za pomocą przycisku `Browse`. Przed uruchomieniem aplikacji, program sprawdza poprawność wejścia m.in. sprawdza czy podano poprawną białkową sekwencję.

W przypadku gdy zostanie podany plik w formacie multiFASTA (np. białko z kilkoma łańcuchami) to algorytm wczyta tylko sekwencję pierwszego łańcucha. To dość często stosowana praktyka w tego typu algorytmach.

Przed uruchomieniem programu jest też sprawdzane czy białko o podanym identyfikatorze (zczytuje go z nagłówka pliku FASTA) znajduje się w bazie DSSP w celu sukcesywnego porównania otrzymanej struktury. Jeżeli nie ma takiej sekwencji użytkownik dostaje powiadomienie o tym. Jeżeli natomiast jest, ale jest zbyt długa względem sekwencji w bazie danych, to program pyta użykownika o zgodę na obcięcie sekwencji do odpowiedniej długości. Gdy nie będzie takiej zgody program uruchomi się, ale statystyki będą mniej wiarygodne.

Następnie należy zaznaczyć pola dla metod, którymi chce się przewidzieć sekwencję.
Kolejna sekcja ustawia plik wynikowy tzn. będą w nim zapisane tylko te informacje, które tutaj zostaną wybrane. Wynik jest zapisywany do pliku z formatem `.output`. Można go potem wczytać bezpośrednio do okna w interfejsie za pomocą przycisku `Open output`.

Zaznaczenie opcji `Scoring` powoduje wypisanie do wynikowego pliku dwóch rzeczy:
1. procent identyczności otrzymanych struktur w stosunku do struktury z bazy DSSP
2. procentowa obecność każdej ze struktur w uzyskanych sekwencjach.

Przycisk `Reset` służy do wyczyszczenia okna wprowadzania.

Aplikację uruchamia się przyciskiem `Run`.
