
# Fait par Alexandre L'Écuyer, 20103530
# Free to use!

def construire_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def rechercher_kmers_dans_base_de_donnees(kmers, base_de_donnees):
    """ Recherche chaque k-mer dans la base de données et renvoie les HSPs initiaux. """
    HSPs = []
    for seq_id, sequence in base_de_donnees.items():
        for kmer in kmers:
            start_pos = sequence.find(kmer)
            while start_pos != -1:
                HSPs.append((seq_id, start_pos, start_pos + len(kmer) - 1, kmer))
                start_pos = sequence.find(kmer, start_pos + 1)
    return HSPs
