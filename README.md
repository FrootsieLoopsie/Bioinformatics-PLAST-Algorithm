
# Bioinformatics-PLAST-Algorithm
BLAST (Basic Local Alignment Search Tool) est une méthode/librairie heuristique de recherche qui permet de trouver très rapidement des régions similaires entre des séquences. Il s'agit d'un des outils les plus utilisés en bioinformatique et biologie comparative. PLAST (Primitive Local Alignment Search Tool) est une version bootleg de cette approche, développée dans le cadre d'études universitaires.

En gros, ce code permet de:
1. Construire tous les kmers présents dans la séquence input pour un k donné.
2. Rechercher les kmers dans les séquences de la base de données en utilisant un algorithme de recherche exact. BLAST utilise des graines qui requièrent des matchs sur chacune des k positions. Les positions où les kmers match représentent les hot spot dénotés HSP (High Scoring Pairs) qui serviront à retrouver l'alignement.
3. Étendre l'alignement, dans les deux directions, aux alentours des HSPs.
4. Évaluer les alignements pour ne garder que ceux pertinents.

C'est particulièrement pratique en bioinformatique pour analyser les séquences génétiques, mais peut avoir des applications pour n'importe quelle comparaison et assemblage de string à partir multiple substrings plus ou moins précis/fiables. J'aimerais appliquer une version optimisée et hardware-accelerated de ce code pour voir à quel point ça pourrait servir dans d'autres sphères que l'academia et les études scientifiques.
