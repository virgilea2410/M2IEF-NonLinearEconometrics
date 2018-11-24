function [] = SetPathDir()
    %Fonction permettant � tout les script Matlab pr�sent dans le m�me dossier
    %que cette fonction SetPathDir(), lors de son appel, d'ajouter les 5
    %chemins des r�pertoires ProjetEconoPART2, Exercices, MyFunctions, Exercice
    %3 et Datas

%Chemin du r�pertoire principal ProjetEconoPART2
DirPath = which('SetPathDir.m');
DirPath = strrep(DirPath, '/SetPathDir.m', '');

%Chemin du r�pertoire Exercices
ExercicesPath = char('/Exercices');
ExercicesPath = strcat(DirPath, ExercicesPath);

%Chemin du r�pertoire MyFunctions
MyFunctionsPath = char('/MyFunctions');
MyFunctionsPath = strcat(ExercicesPath, MyFunctionsPath);

%Chemin du r�pertoire Datas
DatasPath = char('/Datas');
DatasPath = strcat(ExercicesPath, DatasPath);

%Chemin du r�pertoire Exercices 3
Exo3Path = char('/Exercice 3');
Exo3Path = strcat(ExercicesPath, Exo3Path);

%Ajout des chemins dans le compilateur 
addpath(ExercicesPath);
addpath(MyFunctionsPath);
addpath(DatasPath);
addpath(Exo3Path);
addpath(DirPath);

end
