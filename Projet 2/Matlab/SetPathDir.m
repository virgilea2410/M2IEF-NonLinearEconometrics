function [] = SetPathDir()
    %Fonction permettant à tout les script Matlab présent dans le même dossier
    %que cette fonction SetPathDir(), lors de son appel, d'ajouter les 5
    %chemins des répertoires ProjetEconoPART2, Exercices, MyFunctions, Exercice
    %3 et Datas

%Chemin du répertoire principal ProjetEconoPART2
DirPath = which('SetPathDir.m');
DirPath = strrep(DirPath, '/SetPathDir.m', '');

%Chemin du répertoire Exercices
ExercicesPath = char('/Exercices');
ExercicesPath = strcat(DirPath, ExercicesPath);

%Chemin du répertoire MyFunctions
MyFunctionsPath = char('/MyFunctions');
MyFunctionsPath = strcat(ExercicesPath, MyFunctionsPath);

%Chemin du répertoire Datas
DatasPath = char('/Datas');
DatasPath = strcat(ExercicesPath, DatasPath);

%Chemin du répertoire Exercices 3
Exo3Path = char('/Exercice 3');
Exo3Path = strcat(ExercicesPath, Exo3Path);

%Ajout des chemins dans le compilateur 
addpath(ExercicesPath);
addpath(MyFunctionsPath);
addpath(DatasPath);
addpath(Exo3Path);
addpath(DirPath);

end
