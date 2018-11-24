function [] = SetPathExercices()
    %Fonction permettant à tout les script Matlab présent dans le même dossier
    %que cette fonction SetPathDir(), lors de son appel, d'ajouter les 5
    %chemins des répertoires Exercices, MyFunctions, Exercice
    %3 et Datas

%Chemin du répertoire principal Exercices
DirPath = which('SetPathExercices.m');
DirPath = strrep(DirPath, '/SetPathExercices.m', '');

%Chemin du répertoire MyFunctions
MyFunctionsPath = char('/MyFunctions');
MyFunctionsPath = strcat(DirPath, MyFunctionsPath);

%Chemin du répertoire Datas
DatasPath = char('/Datas');
DatasPath = strcat(DirPath, DatasPath);

%Chemin du répertoire Exercice 3
Exo3Path = char('/Exercice 3');
Exo3Path = strcat(DirPath, Exo3Path);

%Ajout des chemins dans le compilateur 
addpath(MyFunctionsPath);
addpath(DatasPath);
addpath(Exo3Path);
addpath(DirPath);

end
