#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#define MAXN 100

static int N, M;
static int adj[MAXN][MAXN];

// Reținem cea mai bună clică globală (și dimensiunea ei)
static int bestClique[MAXN];
static int bestSize = 0;

/********************************************************
 *                 Citire Graf                          *
 ********************************************************/
void readGraph() {
    (void)scanf("%d %d", &N, &M);  // ignorăm return value -> eliminăm warning
    memset(adj, 0, sizeof(adj));
    for (int i = 0; i < M; i++) {
        int x, y;
        (void)scanf("%d %d", &x, &y);  // ignorăm return value -> eliminăm warning
        // presupunem noduri 0..N-1
        adj[x][y] = 1;
        adj[y][x] = 1;
    }
}

/********************************************************
 *              Verificare dacă e clică                 *
 ********************************************************/
int isClique(const int *nodes, int k) {
    for (int i = 0; i < k; i++) {
        for (int j = i+1; j < k; j++) {
            if (!adj[nodes[i]][nodes[j]]) {
                return 0;
            }
        }
    }
    return 1;
}

/********************************************************
 * 1) Backtracking Simplu (Naiv)                        *
 ********************************************************/
void backtrackingNaive(int start, int *current, int sizeCurrent) {
    // Actualizăm best dacă e cazul
    if (sizeCurrent > bestSize) {
        bestSize = sizeCurrent;
        memcpy(bestClique, current, bestSize * sizeof(int));
    }
    // Explorăm în continuare
    for (int v = start; v < N; v++) {
        current[sizeCurrent] = v;
        // Testăm dacă subsetul rămâne clică
        if (isClique(current, sizeCurrent + 1)) {
            backtrackingNaive(v + 1, current, sizeCurrent + 1);
        }
    }
}

/********************************************************
 * 2) Algoritmul Bron–Kerbosch                          *
 ********************************************************/
void bronKerbosch(int *R, int sizeR, int *P, int sizeP, int *X, int sizeX) {
    // Dacă P și X sunt goale, R este o clică maximală
    if (sizeP == 0 && sizeX == 0) {
        if (sizeR > bestSize) {
            bestSize = sizeR;
            memcpy(bestClique, R, bestSize * sizeof(int));
        }
        return;
    }
    
    // Copiem P ca să iterăm în siguranță
    int *Pcopy = (int*)malloc(sizeP * sizeof(int));
    memcpy(Pcopy, P, sizeP * sizeof(int));
    
    for (int i = 0; i < sizeP; i++) {
        int v = Pcopy[i];
        
        // R ∪ {v}
        R[sizeR] = v;
        
        // Construim P \ N(v) și X \ N(v)
        int newP[MAXN], newX[MAXN];
        int sizeNewP = 0, sizeNewX = 0;
        
        for (int j = 0; j < sizeP; j++) {
            int w = P[j];
            if (adj[v][w]) {
                newP[sizeNewP++] = w;
            }
        }
        for (int j = 0; j < sizeX; j++) {
            int w = X[j];
            if (adj[v][w]) {
                newX[sizeNewX++] = w;
            }
        }
        
        // Apel recursiv
        bronKerbosch(R, sizeR + 1, newP, sizeNewP, newX, sizeNewX);
        
        // Mutăm v din P în X
        int pos = -1;
        for (int j = 0; j < sizeP; j++) {
            if (P[j] == v) {
                pos = j;
                break;
            }
        }
        if (pos != -1) {
            P[pos] = P[sizeP - 1];
            sizeP--;
        }
        X[sizeX] = v;
        sizeX++;
    }
    free(Pcopy);
}

/********************************************************
 * 3) Hill Climbing (Local Search)                      *
 ********************************************************/
void hillClimbing(int iterations) {
    srand((unsigned)time(NULL));
    
    int currentClique[MAXN];
    int currentSize = 0;
    
    // Inițializare random
    int startNode = rand() % N;
    currentClique[currentSize++] = startNode;
    
    // Adăugăm nodurile compatibile cu startNode
    for (int v = 0; v < N; v++) {
        if (v == startNode) continue;
        int ok = 1;
        for (int i = 0; i < currentSize; i++) {
            if (!adj[v][currentClique[i]]) {
                ok = 0;
                break;
            }
        }
        if (ok) {
            currentClique[currentSize++] = v;
        }
    }
    
    for (int it = 0; it < iterations; it++) {
        int improved = 0;
        // Încearcă să adaugi un nod care nu e în clică
        for (int v = 0; v < N; v++) {
            // verificăm dacă v e deja în clică
            int inClique = 0;
            for (int i = 0; i < currentSize; i++) {
                if (currentClique[i] == v) {
                    inClique = 1;
                    break;
                }
            }
            if (inClique) continue;
            
            // verificăm dacă v e adiacent la toți din currentClique
            int ok = 1;
            for (int i = 0; i < currentSize; i++) {
                if (!adj[v][currentClique[i]]) {
                    ok = 0;
                    break;
                }
            }
            if (ok) {
                // îl adăugăm
                currentClique[currentSize++] = v;
                improved = 1;
                break; // reluăm din nou
            }
        }
        
        // Dacă nu am găsit nimic de adăugat, încercăm eliminarea unui nod
        if (!improved && currentSize > 1) {
            int removeIndex = rand() % currentSize;
            currentClique[removeIndex] = currentClique[currentSize - 1];
            currentSize--;

            // După eliminare, încercăm 1-2 adăugări noi
            int tries = 2;
            while (tries--) {
                for (int v = 0; v < N; v++) {
                    int inClique = 0;
                    for (int i = 0; i < currentSize; i++) {
                        if (currentClique[i] == v) { 
                            inClique = 1; 
                            break;
                        }
                    }
                    if (inClique) continue;
                    // Compatibil?
                    int ok = 1;
                    for (int i = 0; i < currentSize; i++) {
                        if (!adj[v][currentClique[i]]) {
                            ok = 0;
                            break;
                        }
                    }
                    if (ok) {
                        currentClique[currentSize++] = v;
                        break; 
                    }
                }
            }
        }
        
        // Updatăm best global
        if (currentSize > bestSize) {
            bestSize = currentSize;
            memcpy(bestClique, currentClique, bestSize*sizeof(int));
        }
    }
}

/********************************************************
 * 4) Simulated Annealing (euristică)                   *
 ********************************************************/
void simulatedAnnealing(int maxIter, double startTemp, double alpha) {
    srand((unsigned)time(NULL));
    
    // 1) Inițializare random
    int currentClique[MAXN];
    int currentSize = 0;
    
    // Alegem un nod la întâmplare
    int startNode = rand() % N;
    currentClique[currentSize++] = startNode;
    
    // Adăugăm tot ce e compatibil cu startNode
    for (int v = 0; v < N; v++) {
        if (v == startNode) continue;
        int ok = 1;
        for (int i = 0; i < currentSize; i++) {
            if (!adj[v][currentClique[i]]) {
                ok = 0;
                break;
            }
        }
        if (ok) {
            currentClique[currentSize++] = v;
        }
    }
    
    // Reținem cea mai bună soluție locală
    int bestLocalClique[MAXN];
    int bestLocalSize = currentSize;
    memcpy(bestLocalClique, currentClique, bestLocalSize * sizeof(int));

    // Updatăm best global, dacă e cazul
    if (bestLocalSize > bestSize) {
        bestSize = bestLocalSize;
        memcpy(bestClique, bestLocalClique, bestSize * sizeof(int));
    }
    
    double T = startTemp;
    
    for (int iter = 0; iter < maxIter; iter++) {
        // Construim un vecin (neighbor)
        int neighbor[MAXN];
        memcpy(neighbor, currentClique, currentSize*sizeof(int));
        int neighborSize = currentSize;
        
        // Mișcare: alegem să scoatem un nod sau să adăugăm un nod compatibil
        double r = (double)rand() / RAND_MAX;
        if (r < 0.5 && neighborSize > 0) {
            // scoatem un nod (random)
            int removeIndex = rand() % neighborSize;
            neighbor[removeIndex] = neighbor[neighborSize - 1];
            neighborSize--;
        } else {
            // încercăm să adăugăm un nod care nu e deja în subset
            int candidate = rand() % N;
            // verificăm dacă e deja în subset
            int inClique = 0;
            for (int i = 0; i < neighborSize; i++) {
                if (neighbor[i] == candidate) {
                    inClique = 1;
                    break;
                }
            }
            if (!inClique) {
                // verificăm compatibilitatea cu nodurile din neighbor
                int ok = 1;
                for (int i = 0; i < neighborSize; i++) {
                    if (!adj[candidate][neighbor[i]]) {
                        ok = 0;
                        break;
                    }
                }
                if (ok) {
                    neighbor[neighborSize++] = candidate;
                }
            }
        }
        
        // Calculăm "calitatea" = neighborSize - currentSize
        int delta = neighborSize - currentSize;  
        
        if (delta > 0) {
            // soluția e mai bună, o acceptăm
            memcpy(currentClique, neighbor, neighborSize*sizeof(int));
            currentSize = neighborSize;
        } else {
            // acceptăm cu probabilitate e^(delta/T), (delta<0 => e^(negativ))
            double prob = exp((double)delta / T);
            double rr = (double)rand() / RAND_MAX;
            if (rr < prob) {
                // acceptăm soluția mai slabă
                memcpy(currentClique, neighbor, neighborSize*sizeof(int));
                currentSize = neighborSize;
            }
        }
        
        // Updatăm best local/global
        if (currentSize > bestLocalSize) {
            bestLocalSize = currentSize;
            memcpy(bestLocalClique, currentClique, bestLocalSize*sizeof(int));
            if (bestLocalSize > bestSize) {
                bestSize = bestLocalSize;
                memcpy(bestClique, bestLocalClique, bestSize*sizeof(int));
            }
        }
        
        // Răcim temperatura
        T *= alpha;
        if (T < 1e-8) {
            T = 1e-8;  // evităm zero
        }
    }

    // --- PAS FINAL DE ÎMBUNĂTĂȚIRE LOCALĂ ---
    //  După ce SA se termină, luăm cea mai bună soluție locală
    //  și încercăm să adăugăm noduri compatibile.
    int improved = 1;
    while (improved) {
        improved = 0;
        // căutăm un nod care nu e în bestLocalClique dar se potrivește
        for (int v = 0; v < N; v++) {
            // verificăm dacă v e deja în bestLocalClique
            int inClique = 0;
            for (int i = 0; i < bestLocalSize; i++) {
                if (bestLocalClique[i] == v) {
                    inClique = 1;
                    break;
                }
            }
            if (inClique) continue;

            // verificăm compatibilitatea
            int ok = 1;
            for (int i = 0; i < bestLocalSize; i++) {
                if (!adj[v][bestLocalClique[i]]) {
                    ok = 0;
                    break;
                }
            }
            if (ok) {
                bestLocalClique[bestLocalSize++] = v;
                improved = 1;
                break; // reluăm
            }
        }
    }
    
    // Dacă post-processing a crescut dimensiunea clicei, updatăm best global
    if (bestLocalSize > bestSize) {
        bestSize = bestLocalSize;
        memcpy(bestClique, bestLocalClique, bestSize*sizeof(int));
    }
}

/********************************************************
 * 5) TABU SEARCH (euristică suplimentară, scrisă in C) *
 ********************************************************/
typedef struct {
    int node;
    int isAdd;         // 1 = add, 0 = remove
    int expireIter;    // iteration at which this move expires
} TabuMove;

#define TABU_SIZE 500
#define INF 1000000000

int isTabu(TabuMove *tabuList, int tabuCount, int node, int isAdd, int iter) {
    for (int i = 0; i < tabuCount; i++) {
        if (tabuList[i].node == node &&
            tabuList[i].isAdd == isAdd &&
            tabuList[i].expireIter > iter) 
        {
            return 1; // still tabu
        }
    }
    return 0;
}

void addTabu(TabuMove *tabuList, int *tabuCount,
             int node, int isAdd, int iter, int tabuTenure)
{
    if (*tabuCount < TABU_SIZE) {
        tabuList[*tabuCount].node       = node;
        tabuList[*tabuCount].isAdd      = isAdd;
        tabuList[*tabuCount].expireIter = iter + tabuTenure;
        (*tabuCount)++;
    } else {
        int pos = rand() % TABU_SIZE;
        tabuList[pos].node       = node;
        tabuList[pos].isAdd      = isAdd;
        tabuList[pos].expireIter = iter + tabuTenure;
    }
}

void tabuSearch(int maxIter, int tabuTenure) {
    srand((unsigned)time(NULL));
    
    int currentClique[MAXN];
    int currentSize = 0;
    
    int startNode = rand() % N;
    currentClique[currentSize++] = startNode;
    
    // Add all compatible nodes
    for (int v = 0; v < N; v++) {
        if (v == startNode) continue;
        int ok = 1;
        for (int i = 0; i < currentSize; i++) {
            if (!adj[v][currentClique[i]]) {
                ok = 0;
                break;
            }
        }
        if (ok) {
            currentClique[currentSize++] = v;
        }
    }
    
    // best local
    int bestLocalClique[MAXN];
    int bestLocalSize = currentSize;
    memcpy(bestLocalClique, currentClique, currentSize * sizeof(int));
    
    // update best global
    if (bestLocalSize > bestSize) {
        bestSize = bestLocalSize;
        memcpy(bestClique, bestLocalClique, bestSize * sizeof(int));
    }
    
    // Tabu list
    TabuMove tabuList[TABU_SIZE];
    int tabuCount = 0;
    memset(tabuList, 0, sizeof(tabuList));
    
    for (int iter = 0; iter < maxIter; iter++) {
        int bestNeighborSize = -INF;
        int bestNeighborClique[MAXN];
        
        // Move plan
        int moveNode = -1;
        int moveIsAdd = -1;
        
        // 1) Try ADD moves
        for (int node = 0; node < N; node++) {
            // if already in clique => skip
            int inClique = 0;
            for (int i = 0; i < currentSize; i++) {
                if (currentClique[i] == node) {
                    inClique = 1; 
                    break;
                }
            }
            if (inClique) continue;
            
            // check compatibility
            int ok = 1;
            for (int i = 0; i < currentSize; i++) {
                if (!adj[node][currentClique[i]]) {
                    ok = 0;
                    break;
                }
            }
            if (!ok) continue;
            
            // check Tabu
            if (isTabu(tabuList, tabuCount, node, 1, iter)) {
                continue;
            }
            
            // build neighbor
            int neighborClique[MAXN];
            memcpy(neighborClique, currentClique, currentSize * sizeof(int));
            neighborClique[currentSize] = node;
            int neighborSize = currentSize + 1;
            
            if (neighborSize > bestNeighborSize) {
                bestNeighborSize = neighborSize;
                memcpy(bestNeighborClique, neighborClique, neighborSize*sizeof(int));
                moveNode  = node;
                moveIsAdd = 1; // add
            }
        }
        
        // 2) Try REMOVE moves
        for (int idx = 0; idx < currentSize; idx++) {
            int node = currentClique[idx];
            // check Tabu
            if (isTabu(tabuList, tabuCount, node, 0, iter)) {
                continue;
            }
            // build neighbor
            int neighborClique[MAXN];
            memcpy(neighborClique, currentClique, currentSize*sizeof(int));
            neighborClique[idx] = neighborClique[currentSize - 1];
            int neighborSize = currentSize - 1;
            
            if (neighborSize > bestNeighborSize) {
                bestNeighborSize = neighborSize;
                memcpy(bestNeighborClique, neighborClique, neighborSize*sizeof(int));
                moveNode  = node;
                moveIsAdd = 0; // remove
            }
        }
        
        // no valid neighbor?
        if (bestNeighborSize == -INF) {
            break;
        }
        
        // apply best neighbor
        memcpy(currentClique, bestNeighborClique, bestNeighborSize*sizeof(int));
        currentSize = bestNeighborSize;
        
        // add move to Tabu list
        addTabu(tabuList, &tabuCount, moveNode, moveIsAdd, iter, tabuTenure);
        
        // update best local/global
        if (currentSize > bestLocalSize) {
            bestLocalSize = currentSize;
            memcpy(bestLocalClique, currentClique, bestLocalSize*sizeof(int));
            if (bestLocalSize > bestSize) {
                bestSize = bestLocalSize;
                memcpy(bestClique, bestLocalClique, bestSize * sizeof(int));
            }
        }
    }
}

/********************************************************
 *        Funcție de afișare a soluției curente         *
 ********************************************************/
void printSolution(double elapsedTime, int runIndex) {
    // Afișăm: indexul rularii, timpul, dimensiunea clicei, nodurile
    printf("Rulare %3d | Timp: %.6f s | Dimensiune: %d | Noduri: ",
           runIndex, elapsedTime, bestSize);
    for (int i = 0; i < bestSize; i++) {
        printf("%d ", bestClique[i]);
    }
    printf("\n");
    // Mică pauză între afișări
    // sleep(1); // optional
}

/********************************************************
 *                       main()                          *
 ********************************************************/
int main() {
    readGraph();
    
    // Vom rula fiecare dintre cele 5 metode de 500 de ori
    // => 5 * 500 = 1000 rânduri de output.
    // In plus, la finalul celor 500 de rulari de fiecare tip,
    // afisam timpul mediu si nr. de raspunsuri pt. fiecare dimensiune.

    // Bufor pentru submulțimea temporară (în backtracking)
    int *temp = (int*)malloc(N * sizeof(int));
    
    // Pentru stocarea timpilor si a dimensiunilor
    static double timesBack[500];
    static int    sizesBack[500];
    
    static double timesBron[500];
    static int    sizesBron[500];
    
    static double timesHill[500];
    static int    sizesHill[500];
    
    static double timesSA[500];
    static int    sizesSA[500];
    
    static double timesTabu[500];
    static int    sizesTabu[500];
    
    // 1) BACKTRACKING SIMPLU
    printf("=== BACKTRACKING SIMPLU ===\n");
    for(int run = 0; run < 500; run++) {
        clock_t start = clock();
        
        bestSize = 0;  // resetăm înainte de fiecare rulare
        backtrackingNaive(0, temp, 0);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        // Afișăm pe ecran această rulare (opțional)
        printSolution(elapsed, run+1);
        
        // Salvăm în array-urile noastre
        timesBack[run] = elapsed;
        sizesBack[run] = bestSize;
        
        // Reset bestClique intre rulări
        for(int i = 0; i < bestSize; i++) {
            bestClique[i] = 0; 
        }
        bestSize = 0;  
    }
    // La final, calculăm și afișăm timpul mediu + frecvența dimensiunilor
    {
        double sum = 0.0;
        for (int i = 0; i < 500; i++) {
            sum += timesBack[i];
        }
        double avgTime = sum / 500.0;
        
        // frecvența pentru fiecare dimensiune
        int freq[101];
        memset(freq, 0, sizeof(freq));
        // (presupunem că nicio clică nu va depăși 100, deoarece N<=100)
        for (int i = 0; i < 500; i++) {
            freq[sizesBack[i]]++;
        }
        
        printf("=== Rezultate finale Backtracking ===\n");
        printf("Timp mediu: %.6f s\n", avgTime);
        printf("Dimensiune --> NumarRulari\n");
        for (int d = 0; d <= 100; d++) {
            if (freq[d] > 0) {
                printf("%d --> %d\n", d, freq[d]);
            }
        }
    }
    
    // 2) BRON–KERBOSCH
    printf("\n=== BRON–KERBOSCH ===\n");
    for(int run = 0; run < 500; run++) {
        clock_t start = clock();
        
        bestSize = 0;
        int R[MAXN], P[MAXN], X[MAXN];
        for(int i = 0; i < N; i++) {
            P[i] = i;
        }
        bronKerbosch(R, 0, P, N, X, 0);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printSolution(elapsed, run+1);
        
        timesBron[run] = elapsed;
        sizesBron[run] = bestSize;
        
        for(int i = 0; i < bestSize; i++) {
            bestClique[i] = 0; 
        }
        bestSize = 0;  
    }
    // Rezultate finale Bron–Kerbosch
    {
        double sum = 0.0;
        for (int i = 0; i < 500; i++) {
            sum += timesBron[i];
        }
        double avgTime = sum / 500.0;
        int freq[101];
        memset(freq, 0, sizeof(freq));
        for (int i = 0; i < 500; i++) {
            freq[sizesBron[i]]++;
        }
        
        printf("=== Rezultate finale Bron-Kerbosch ===\n");
        printf("Timp mediu: %.6f s\n", avgTime);
        printf("Dimensiune --> NumarRulari\n");
        for (int d = 0; d <= 100; d++) {
            if (freq[d] > 0) {
                printf("%d --> %d\n", d, freq[d]);
            }
        }
    }
    
    // 3) HILL CLIMBING
    printf("\n=== HILL CLIMBING ===\n");
    for(int run = 0; run < 500; run++) {
        clock_t start = clock();
        
        bestSize = 0;
        hillClimbing(500); // ex: 500 iterații
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printSolution(elapsed, run+1);
        
        timesHill[run] = elapsed;
        sizesHill[run] = bestSize;
        
        for(int i = 0; i < bestSize; i++) {
            bestClique[i] = 0; 
        }
        bestSize = 0;  
        usleep(800000);
    }
    // Rezultate finale Hill Climbing
    {
        double sum = 0.0;
        for (int i = 0; i < 500; i++) {
            sum += timesHill[i];
        }
        double avgTime = sum / 500.0;
        int freq[101];
        memset(freq, 0, sizeof(freq));
        for (int i = 0; i < 500; i++) {
            freq[sizesHill[i]]++;
        }
        
        printf("=== Rezultate finale Hill Climbing ===\n");
        printf("Timp mediu: %.6f s\n", avgTime);
        printf("Dimensiune --> NumarRulari\n");
        for (int d = 0; d <= 100; d++) {
            if (freq[d] > 0) {
                printf("%d --> %d\n", d, freq[d]);
            }
        }
    }
    
    // 4) SIMULATED ANNEALING
    printf("\n=== SIMULATED ANNEALING ===\n");
    for(int run = 0; run < 500; run++) {
        clock_t start = clock();
        
        bestSize = 0;
        simulatedAnnealing(/*maxIter*/ 1500, /*startTemp*/ 10.0, /*alpha*/ 0.98);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printSolution(elapsed, run+1);
        
        timesSA[run] = elapsed;
        sizesSA[run] = bestSize;
        
        for(int i = 0; i < bestSize; i++) {
            bestClique[i] = 0; 
        }
        bestSize = 0; 
        usleep(800000); 
    }
    // Rezultate finale SA
    {
        double sum = 0.0;
        for (int i = 0; i < 500; i++) {
            sum += timesSA[i];
        }
        double avgTime = sum / 500.0;
        int freq[101];
        memset(freq, 0, sizeof(freq));
        for (int i = 0; i < 500; i++) {
            freq[sizesSA[i]]++;
        }
        
        printf("=== Rezultate finale Simulated Annealing ===\n");
        printf("Timp mediu: %.6f s\n", avgTime);
        printf("Dimensiune --> NumarRulari\n");
        for (int d = 0; d <= 100; d++) {
            if (freq[d] > 0) {
                printf("%d --> %d\n", d, freq[d]);
            }
        }
    }

    // 5) TABU SEARCH
    printf("\n=== TABU SEARCH ===\n");
    for(int run = 0; run < 500; run++) {
        clock_t start = clock();
        
        bestSize = 0;
        // param ex: maxIter=1000, tabuTenure=10
        tabuSearch(1000, 10);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printSolution(elapsed, run+1);
        
        timesTabu[run] = elapsed;
        sizesTabu[run] = bestSize;
        
        for(int i = 0; i < bestSize; i++) {
            bestClique[i] = 0; 
        }
        bestSize = 0; 
        usleep(800000); 
    }
    // Rezultate finale Tabu
    {
        double sum = 0.0;
        for (int i = 0; i < 500; i++) {
            sum += timesTabu[i];
        }
        double avgTime = sum / 500.0;
        int freq[101];
        memset(freq, 0, sizeof(freq));
        for (int i = 0; i < 500; i++) {
            freq[sizesTabu[i]]++;
        }
        
        printf("=== Rezultate finale Tabu Search ===\n");
        printf("Timp mediu: %.6f s\n", avgTime);
        printf("Dimensiune --> NumarRulari\n");
        for (int d = 0; d <= 100; d++) {
            if (freq[d] > 0) {
                printf("%d --> %d\n", d, freq[d]);
            }
        }
    }
    
    free(temp);
    return 0;
}
