from Bio import SeqIO
import math
import argparse


# Par Alexandre L'Écuyer 20103530

# Comment rouler ce fichier:
#   1 - Installer les librairies, avec: pip install biopython
#   2 - Placer le fichier fasta au même répertoire que ce fichier.
#   3 - Exécuter ce fichier avec python, avec les arguments valides.


# ---------------------------------------------------------------------------------------
# LECTURE DU FICHIER
def read_fasta(chemin_fichier):
    """ Utilise SeqIO de Biopython pour lire un fichier FASTA et construire la base de données/dictionnaire. """
    database = {}
    for read in SeqIO.parse(chemin_fichier, "fasta"):
        database[read.id] = str(read.seq)
    return database


# ---------------------------------------------------------------------------------------
# SEEDS AND KMER DETECTION:
# Pour chacune des séquences de la banque de données, on va trouver les positions des occurences
# de chaque kmer (substring de taille k). Ces positions (et l'alignement correspondant) représentent
# les HSPs initiaux qu'il faudrait ensuite étendre dans les deux sens.
def construire_kmers(sequence, k):
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def recherche_exacte_kmers(kmers, database):
    """ Recherche chaque k-mer dans la base de données et renvoie les HSPs initiaux avec les positions dans sequence. """
    HSPs = []
    for seq_id, db_seq in database.items():
        for kmer in kmers:
            start_pos_db = db_seq.find(kmer)
            while start_pos_db != -1:
                # Trouver la position correspondante dans 'sequence'
                end_pos_db = start_pos_db + len(kmer) - 1
                HSPs.append((seq_id, start_pos_db, end_pos_db, kmer))
                start_pos_db = db_seq.find(kmer, start_pos_db + 1)
    return HSPs


# ---------------------------------------------------------------------------------------
# GREEDY HSP EXTENSION HEURISTIC
# Il faut ensuite étendre la similitude dans les deux sens (en ignorant la possibilité de gaps) le long
# de chaque séquence, à partir du HSP initial, de manière à ce que le score cumulé puisse être amélioré.
# Seule l'extension produisant le meilleur score pour le HSP étendu sera choisi.
def extend_hsp(sequence, hsp, database, score_match=5, score_mismatch=-4, seuil_extension=4):
    seq_id, db_left, db_right, kmer = hsp
    db_seq = database[seq_id]

    # Initialisez le meilleur score au score de l'overlap initial.
    score_total = score_max = (db_right - db_left + 1) * score_match
    sq_left = sequence.find(kmer)
    sq_right = sq_left + len(kmer) - 1
    best_extension = (seq_id, db_left, db_right, sq_left, sq_right, score_max)

    while True:
        score_left, score_right = 0, 0

        # Vérifiez l'extension à gauche
        if (db_left > 0) and (sq_left > 0):
            score_left = score_match if sequence[sq_left - 1] == db_seq[db_left - 1] else score_mismatch

        # Vérifiez l'extension à droite
        if (db_right < len(db_seq) - 1) and (sq_right < len(sequence) + 1):
            score_right = score_match if sequence[sq_right + 1] == db_seq[db_right + 1] else score_mismatch

        # Arrêter si aucune extension n'est possible
        if score_left == 0 and score_right == 0:
            break

        # Choisir la meilleure extension
        if (score_left != 0) and (score_left >= score_right):
            db_left -= 1
            sq_left -= 1
            score_total += score_left

        elif (score_right != 0):
            db_right += 1
            sq_right += 1
            score_total += score_right

        else:
            return best_extension


        # Mettre à jour le meilleur score
        if score_total > score_max:
            score_max = score_total
            best_extension = (seq_id, db_left, db_right, sq_left, sq_right, score_max)

        # Si le score total baisse trop, revenir à la meilleure extension et arrêter
        elif score_max - score_total >= seuil_extension:
            break

    return best_extension



# ---------------------------------------------------------------------------------------
# MERGE ALIGNED HSP MATCHES
# Fusionner tous les HSPs se chevauchant et impliquant les mêmes séquences
def fusionner_hsp(hsps, database, sequence):
    """ Fusionne les HSPs chevauchants. """
    if(len(hsps) == 0):
        return None

    # Tri des HSPs par identifiant de séquence et position de début
    hsps.sort(key=lambda x: (x[0], x[1]))
    fusionnes = []

    # Initialiser le premier HSP comme étant le HSP courant à comparer
    hsp_courant = hsps[0]
    for hsp in hsps[1:]:
        # Vérifier si le HSP actuel et le HSP courant se chevauchent et concernent la même séquence
        if hsp[0] == hsp_courant[0] and hsp[1] <= hsp_courant[2]:
            # Si les HSPs se chevauchent, mettre à jour la fin du HSP courant si nécessaire
            hsp_courant = (hsp_courant[0], hsp_courant[1], max(hsp_courant[2], hsp[2]), hsp_courant[3], max(hsp_courant[4], hsp[4]), hsp_courant[5])
        else:
            # Si les HSPs ne se chevauchent pas, ajouter le HSP courant à la liste des fusionnés, puis passer au prochain
            fusionnes.append(hsp_courant)
            hsp_courant = hsp

    # Ajouter le dernier HSP courant à la liste des fusionnés
    fusionnes.append(hsp_courant)
    return fusionnes


# ---------------------------------------------------------------------------------------
# HSP STATISTICS
# Pour comparer puis filtrer HSPs selon un système de scorage robuste. Plutôt que comparer en fct de leurs scores bruts,
# puisque ces derniers sont biaisés par la taille des alignements, j'utilise un bitscore établi dans le cours de BioInfo IFT.
# À partir de ce score vous pouvez calculer la e-value qui représente le nombre attendu, par chance, de HSP avec un bitscore d'au moins B.
# Plus le bitscore est grand (et plus la e-value est petite), plus le HSP est pertinent -> haut score de signifiance.

def calculer_bitscore(score_brut, lambda_value, K):
    """ Calcule le bitscore à partir du score brut. """
    return round((lambda_value * score_brut - math.log(K)) / math.log(2))


def calculer_e(m, n, bitscore):
    """ Calcule la e-value pour un bitscore donné. """
    return m * n * 2 ** (-bitscore)


def filtrer_hsp_significatifs(hsps, size_totale_reads, size_of_sequence_to_lookup, seuil):
    """ Filtre les HSPs en fonction de leur e-value. """
    hsp_significatifs = []
    if(hsps != None and len(hsps) > 0):
        for hsp in hsps:
            bitscore = calculer_bitscore(hsp[5], lambda_value=0.192, K=0.176)
            e_value = calculer_e(size_totale_reads, size_of_sequence_to_lookup, bitscore)

            if e_value < seuil:
                hsp_significatifs.append(hsp + (e_value,))
    return hsp_significatifs


if __name__ == "__main__":
    # Définir l'analyseur d'arguments de la ligne de commande
    parser = argparse.ArgumentParser(description='PLAST: Primitive Local Alignment Search Tool.')
    parser.add_argument('-i', type=str, required=True, help='Input nucleotide sequence to search.')
    parser.add_argument('-db', type=str, required=True, help='Database file in FASTA format.')
    parser.add_argument('-E', type=int, default=4, help='Extension threshold (default: 4).')
    parser.add_argument('-ss', type=float, default=1e-3, help='Significance score threshold (default: 1e-3).')
    parser.add_argument('-seed', type=str, default='11111111111', help='Seed pattern for k-mers (default: "11111111111").')

    # Parser les arguments
    args = parser.parse_args()
    sequence_to_search = args.i
    chemin_fichier = args.db
    seuil_extension = args.E
    seuil_e_value = args.ss
    seed_pattern = args.seed

    # La taille des k-mers est la longueur du motif de la graine
    k_size = len(seed_pattern)  # Vous pouvez ajuster la taille des k-mers ici

    database = read_fasta(chemin_fichier)
    kmers = construire_kmers(sequence_to_search, k_size)
    HSPs_initiaux = recherche_exacte_kmers(kmers, database)
    HSPs_etendus = [extend_hsp(sequence_to_search, hsp, database, score_match=5, score_mismatch=-4, seuil_extension=seuil_extension)  for hsp in HSPs_initiaux]
    HSPs_fusionnes = fusionner_hsp(HSPs_etendus, database, sequence_to_search)

    # Trier les HSPs par e-value croissante
    HSPs_significatifs = filtrer_hsp_significatifs(HSPs_fusionnes, len(database), len(sequence_to_search), seuil_e_value)
    HSPs_significatifs.sort(key=lambda x: x[-1])

    # Sélectionner le meilleur HSP pour chaque séquence
    meilleurs_hsp_par_seq = {}
    for hsp in HSPs_significatifs:
        seq_id = hsp[0]
        if seq_id not in meilleurs_hsp_par_seq or hsp[-1] < meilleurs_hsp_par_seq[seq_id][-1]:
            meilleurs_hsp_par_seq[seq_id] = hsp

    # Affichage des résultats
    for seq_id, hsp in meilleurs_hsp_par_seq.items():
        print(f">{seq_id}")
        print(f"# Best HSP score: {hsp[5]}, bitscore: {calculer_bitscore(hsp[5], 0.192, 0.176)}, evalue: {hsp[-1]:.2e}")

        start, end, start_seq, end_seq = hsp[1], hsp[2], hsp[3], hsp[4]
        print(f"{start} {database[seq_id][start:end + 1]} {end}")
        print(f"{start_seq} {sequence_to_search[start_seq:end_seq + 1]} {end_seq}\n")

    print("--------------------------------------------------")
    print(f"Total : {len(meilleurs_hsp_par_seq)}")
