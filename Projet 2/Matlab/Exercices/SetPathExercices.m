function [] = SetPathExercices()
    %Fonction permettant � tout les script Matlab pr�sent dans le m�me dossier
    %que cette fonction SetPathDir(), lors de son appel, d'ajouter les 5
    %chemins des r�pertoires Exercices, MyFunctions, Exercice
    %3 et Datas

%Chemin du r�pertoire principal Exercices
DirPath = which('SetPathExercices.m');
DirPath = strrep(DirPath, '/SetPathExercices.m', '');

%Chemin du r�pertoire MyFunctions
MyFunctionsPath = char('/MyFunctions');
MyFunctionsPath = strcat(DirPath, MyFunctionsPath);

%Chemin du r�pertoire Datas
DatasPath = char('/Datas');
DatasPath = strcat(DirPath, DatasPath);

%Chemin du r�pertoire Exercice 3
Exo3Path = char('/Exercice 3');
Exo3Path = strcat(DirPath, Exo3Path);

%Ajout des chemins dans le compilateur 
addpath(MyFunctionsPath);
addpath(DatasPath);
addpath(Exo3Path);
addpath(DirPath);

end
