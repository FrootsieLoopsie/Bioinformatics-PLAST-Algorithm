
from Bio import SeqIO

# Par Alexandre L'Écuyer 20103530

PRINT_RESULTS = False

# Comment rouler ce fichier:
#   1 - Installer les librairies, avec:
#           pip install biopython
#           pip install numpy

#   2 - Installer le package graphviz sur l'ordi, si pas déjà fait:
#           Ubuntu: sudo apt-get install graphviz
#           MacOS: brew install graphviz
#           Windows: https://stackoverflow.com/a/44005139 ou https://graphviz.org/download/

#   3 - Placer reads.fq au même répertoire que ce fichier.
#   4 - Exécuter ce fichier avec python.


#---------------------------------------------------------------------------------------
# LECTURE DU FICHIER
def lire_fichier_fasta(chemin_fichier):
    """ Utilise SeqIO de Biopython pour lire un fichier FASTA et construire la base de données. """
    base_de_donnees = {}
    for read in SeqIO.parse(chemin_fichier, "fasta"):
        base_de_donnees[read.id] = str(read.seq)
    return base_de_donnees


#---------------------------------------------------------------------------------------
# SEEDS AND KMER DETECTION:
# Pour chacune des séquences de la banque de données, on va trouver les positions des occurences
# de chaque kmer (substring de taille k). Ces positions (et l'alignement correspondant) représentent
# les HSPs initiaux qu'il faudrait ensuite étendre dans les deux sens.
def construire_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def rechercher__exacte_kmers(kmers, base_de_donnees):
    """ Recherche chaque k-mer dans la base de données et renvoie les HSPs initiaux. """
    HSPs = []
    for seq_id, sequence in base_de_donnees.items():
        for kmer in kmers:
            start_pos = sequence.find(kmer)
            while start_pos != -1:
                HSPs.append((seq_id, start_pos, start_pos + len(kmer) - 1, kmer))
                start_pos = sequence.find(kmer, start_pos + 1)
    return HSPs


#---------------------------------------------------------------------------------------
# GREEDY HSP EXTENSION HEURISTIC
# Il faut ensuite étendre la similitude dans les deux sens (en ignorant la possibilité de gaps) le long
# de chaque séquence, à partir du HSP initial, de manière à ce que le score cumulé puisse être amélioré.
# Seule l'extension produisant le meilleur score pour le HSP étendu sera choisi.

def etendre_hsp(sequence, hsp, base_de_donnees, score_match=5, score_mismatch=-4, seuil_extension=4):
    """ Étend un HSP en maximisant le score d'alignement avec l'heuristique gloutonne. """
    seq_id, start, end = hsp
    db_seq = base_de_donnees[seq_id]
    score_max = 0
    best_extension = (start, end)

    # Initialiser les variables pour l'extension
    score_cumule = 0
    ext_gauche, ext_droite = start, end

    while True:
        score_gauche = score_droite = 0
        # Vérifier l'extension à gauche
        if ext_gauche > 0 and sequence[ext_gauche - 1] == db_seq[ext_gauche - 1]:
            score_gauche = score_match
        else:
            score_gauche = score_mismatch

        # Vérifier l'extension à droite
        if ext_droite < len(db_seq) - 1 and sequence[ext_droite + 1] == db_seq[ext_droite + 1]:
            score_droite = score_match
        else:
            score_droite = score_mismatch

        # Choisir la meilleure extension
        if score_gauche >= score_droite:
            score_cumule += score_gauche
            ext_gauche -= 1
        else:
            score_cumule += score_droite
            ext_droite += 1

        # Mettre à jour le score maximum et la meilleure extension
        if score_cumule > score_max:
            score_max = score_cumule
            best_extension = (ext_gauche, ext_droite)
        elif score_max - score_cumule >= seuil_extension:
            # Arrêter l'extension si le score cumulé descend en dessous du seuil
            break
        elif ext_gauche == 0 or ext_droite == len(db_seq) - 1:
            # Arrêter l'extension si on atteint la fin de l'une des séquences
            break

    # Retourner le HSP étendu avec le meilleur score
    return seq_id, best_extension[0], best_extension[1], score_max