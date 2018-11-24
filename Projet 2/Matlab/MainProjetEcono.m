%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projet Econom�trie Non Lin�aire, Partie 2 : Markov Switching AR Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nettoyage de la Ligne de Commande, du Workspace, fermeture des
%Graphiques, d�sacticvations des Warnings Matlab
clear all;
close all;
clc;
warning('off');

%Ajout des chemins de nos dossiers dans le compilateur pour qu'il puisse
%trouver les fonctions lors de la compilation
SetPathDir();

%Execution de chaque exercice du projet
Exercice1;
Exercice2;
Exercice3;