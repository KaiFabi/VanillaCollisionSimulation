#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <float.h>

/***************/
/* Definitions */
/***************/
#define MAX(a,b) ( (a)>(b)?(a):(b) )
#define MIN(a,b) ( (a)<(b)?(a):(b) )

#define DIM 2
#define EPS_SOFT_LENGTH 0.0001
#define SIGMA 1.0
#define EPSILON 5.0         // binding force
#define RCUT 2.5 * SIGMA
#define PARTICLE_DISTANCE 1.122462048309373 * SIGMA // mesh size

#define HYPOT2(x, y) sqrt(sqr(x) + sqr(y))
#define sqr(x) ((x)*(x))

#define TRUE 1
#define FALSE 0
#define POWDIM 4

typedef enum{particle, pseudoParticle} nodetype;

/**************/
/* Structures */
/**************/
typedef struct{
    double m;
    double x[DIM];
    double v[DIM];
    double F[DIM];
    double F_old[DIM];
    int moved;
    int todelete;
} Particle;

typedef struct{
    double lower[DIM];
    double upper[DIM];
} Box;

typedef struct TreeNode{
    Particle p;
    Box box;
    nodetype node;
    struct TreeNode *son[POWDIM];
} TreeNode;

/**************/
/* Prototypes */
/**************/
void outputResults(TreeNode *root, FILE *fp);
void freeTree(TreeNode *root);
void moveParticles_BH(TreeNode *root);
void setFlags(TreeNode *t);
void moveLeaf(TreeNode *t, TreeNode *root);
void repairTree(TreeNode *t);
void timeIntegration(double t, double dt, double Tmax, TreeNode *root, Box box);
void initData(TreeNode **root, Box *domain, int N);
void inputParameters(double *dt, double *Tmax, Box *box, int *N);
void insertTree(Particle *p, TreeNode *t);
void updateX(Particle *p, double dt);
void compX(TreeNode *root, double dt);
void compV(TreeNode *root, double dt);
void compF(TreeNode *t, TreeNode *root);
void genData(Particle *p, int n);
int outsideBox(TreeNode *t);
int sonNumber(Box *box, Box *sonbox, Particle *p);
double urand(double, double);

/********/
/* Main */
/********/
int main(int argc, char *argv[]){
    TreeNode *root;
    Box box;
    double dt;
    double Tmax;
    int N;
    inputParameters(&dt, &Tmax, &box, &N);
    initData(&root, &box, N);
    timeIntegration(0, dt, Tmax, root, box);
    freeTree(root);
    return 0;
}

/*************/
/* Functions */
/*************/

void outputResults(TreeNode *t, FILE *fp){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            outputResults(t->son[i], fp);
        }
        // Operation on *t
        if(t->node == particle){
            fprintf(fp, "%lf %lf %lf %lf\n",t->p.x[0],t->p.x[1],t->p.v[0],t->p.v[1]);
        }
    }
}

void moveParticles_BH(TreeNode *root){
    setFlags(root);
    moveLeaf(root, root);
    repairTree(root);
}

void setFlags(TreeNode *t){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            setFlags(t->son[i]);
        }
    // Operation on *t
        t->p.moved = FALSE;
        t->p.todelete = FALSE;
    }
}

// Check if particle is still in the box
int outsideBox(TreeNode *t){
    for(int d=0; d<DIM; d++){
        if((t->p.x[d] < t->box.lower[d]) || (t->p.x[d] > t->box.upper[d])){
            return TRUE; // return TRUE if particle is outside of his box
        }
    }
    return FALSE;
}

// Move particels in tree if they are outside of their boxes
void moveLeaf(TreeNode *t, TreeNode *root){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++)
                moveLeaf(t->son[i], root);
        // Operation on *t
        if((t->node == particle)&&(t->p.moved == FALSE)){
            t->p.moved = TRUE;
            if(outsideBox(t) == TRUE){
                insertTree(&t->p, root);
                t->p.todelete = TRUE;
            }
        }
    }
}

// Repair tree
void repairTree(TreeNode *t){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            repairTree(t->son[i]);
        }
        // Operation on *t
        if(t->node != particle){
            int numberofsons = 0;
            int d = -1;
            for(int i=0; i<POWDIM; i++) {
                if(t->son[i] != NULL) {
                    if(t->son[i]->p.todelete == TRUE){
                        free(t->son[i]);
                        t->son[i] = NULL;
                    }
                    else {
                        numberofsons++;
                        d=i;
                    }
                }
            }

            if(numberofsons == 0){
                t->p.todelete = TRUE;
            }else if ((numberofsons == 1) && (t->son[d]->node == particle)){
                t->p = t->son[d]->p;
                t->node = t->son[d]->node;
                free(t->son[d]);
                t->son[d] = NULL;
            }

        }
    }
}

void freeTree(TreeNode *t){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            freeTree(t->son[i]);
        }
        // Operation on *t
        free(t);
    }
}

double distanceDomainBoxParticle2D(TreeNode *td, TreeNode *t){
    double x_min = td->box.lower[0];
    double x_max = td->box.upper[0];
    double y_min = td->box.lower[1];
    double y_max = td->box.upper[1];
    double x = t->p.x[0];
    double y = t->p.x[1];

    if (x < x_min) {
        if (y < y_min)
            return HYPOT2(x_min-x, y_min-y);
        else if (y <= y_max)
            return x_min - x;
        else
            return HYPOT2(x_min-x, y_max-y);
    } else if (x <= x_max) {
        if (y < y_min)
            return y_min - y;
        else if (y <= y_max)
            return 0;
        else
            return y - y_max;
    } else {
        if (y < y_min)
            return HYPOT2(x_max-x, y_min-y);
        else if (y <= y_max)
            return x - x_max;
        else
            return HYPOT2(x_max-x, y_max-y);
    }
}


void updateX(Particle *p, double dt){
    double a = dt*0.5/p->m;
    for(int d=0; d<DIM; d++){
        p->x[d] += dt*(p->v[d] + a*p->F[d]);
        p->F_old[d] = p->F[d];
    }
};

void updateV(Particle *p, double dt){
    double a = dt*0.5/p->m;
    for(int d=0; d<DIM; d++){
        p->v[d] += a*(p->F[d] + p->F_old[d]);
    }
};

void compX(TreeNode *t, double dt){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compX(t->son[i], dt);
        }
        // Operationen on *t
        if(t->node == particle){
            updateX(&t->p, dt);
        }
    }
}

void compV(TreeNode *t, double dt){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compV(t->son[i], dt);
        }
        // Operationen on *t
        if(t->node == particle){
            updateV(&t->p, dt);
        }
    }
}

// Lennard-Jones Force between two Particles
void force(Particle *i, Particle *j){

    double r=0;
    for(int d=0; d<DIM; d++){
        r += sqr(j->x[d] - i->x[d]);
    }
    r += EPS_SOFT_LENGTH;
    double s = sqr(SIGMA)/ r;
    s = sqr(s) * s;
    double f = 24 * EPSILON * s / r * (1 - 2*s);

    for(int d=0; d<DIM; d++){
        i->F[d] += f*(j->x[d] - i->x[d]);
    }
}

void force_tree(TreeNode *tl, TreeNode *t, int level){
    if((t != tl)&&(t != NULL)){
        if( (distanceDomainBoxParticle2D(t, tl) < RCUT) || (level == 0)){
            if(t->node == particle){
                force(&tl->p, &t->p);
            }else{
                for(int i=0; i<POWDIM; i++){
                    force_tree(tl, t->son[i], level+1);
                }
            }
        }
    }
}

void compF(TreeNode *t, TreeNode *root){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compF(t->son[i], root);
        }
        // Operation on *t
        if(t->node == particle){
            for(int d=0; d<DIM; d++){
                t->p.F[d] = 0;
            }
            force_tree(t, root, 0);
        }
    }
}


void get_file_name(int k, char* buffer, size_t buflen){
    snprintf(buffer, buflen, "logs/data_eps5_%d.log", k);
}

void outputResults2File(TreeNode *root, int n){
    const size_t BUFLEN = 50;
    char file_name[BUFLEN];
    get_file_name(n, file_name, BUFLEN);
    FILE *fp;
    fp = fopen(file_name, "w+");
    outputResults(root, fp);
    fclose(fp);
}

void timeIntegration(double t, double dt, double Tmax, TreeNode *root, Box box){
    int n=0;
    while (t<Tmax){
        t += dt;

        if(n%50 == 0){
            outputResults2File(root, n);    // print some results
        }

        compX(root, dt);
        compF(root, root);
        compV(root, dt);
        moveParticles_BH(root);
        n++;
    }
}

void initData(TreeNode **root, Box *domain, int N){
    // Generate Data
    Particle *p = (Particle*)malloc(N*sizeof(*p)); 
    genData(p, N);

    *root = (TreeNode*)calloc(1, sizeof(TreeNode));
    (*root)->p = p[0];
    (*root)->box = *domain;

    for(int i=1; i<N; i++){
        insertTree(&p[i], *root);
    }
    free(p);
}

void inputParameters(double *dt, double *Tmax, Box *box, int *N){
    *dt = 0.00005;
    *Tmax = 10;
    *N = 10100;
    for(int d=0; d<DIM; d++){
        box->lower[d] = -1000;
        box->upper[d] = 1000;
    }
}

double urand(double low, double high){
    return low+((double)rand() / (double)RAND_MAX)*(high-low);
}

void genData(Particle *p, int n){
    time_t t;
    srand((unsigned) time(&t));

    int N_1 = 100;
    int N_2 = 10000;

    double dx = 1.01 * PARTICLE_DISTANCE;
    double dy = 1.01 * PARTICLE_DISTANCE;

    // Object 1
    double px = 0;
    double py = 0;
    double pos_x = -9.*dx;
    double pos_y = 0.5 * (dy * 100);

    int c = 0;
    for(int i=0; i<N_1; i++){
        c++;
        p[i].x[0] = pos_x + px;
        p[i].x[1] = pos_y + py;
        px += dx;
        if(c%10 == 0){
            px = 0;
            py += dy;
        }
        p[i].v[0] = 100 + urand(-0.001, 0.001);
        p[i].v[1] = 0 + urand(-0.001, 0.001);
        p[i].m = 1;
    }

    // Object 2
    px = 0;
    py = 0;
    pos_x = 2 * dx;
    pos_y = 0;
    c = 0;
    int k = 0;
    for(int i=N_1; i<N_1+N_2; i++){
        c++;
        p[i].x[0] = pos_x + px;
        p[i].x[1] = pos_y + py;
        py += dy;
        if(c%100 == 0){
            k++;
            px += dx;
            py = 0;
            if(k%10 == 0){
                px += 5*dx;
            }
        }
        p[i].v[0] = urand(-0.001, 0.001); 
        p[i].v[1] = urand(-0.001, 0.001); 
        p[i].m = 1;
    }
}

int sonNumber(Box *box, Box *sonbox, Particle *p){
    int b=0;
    for(int d=DIM-1; d>=0; d--){
        if(p->x[d] < .5*(box->upper[d] + box->lower[d])){
            b = 2*b;
            sonbox->lower[d] = box->lower[d];
            sonbox->upper[d] = .5*(box->upper[d] + box->lower[d]);
        }else{
            b = 2*b+1;
            sonbox->lower[d] = .5*(box->upper[d] + box->lower[d]);
            sonbox->upper[d] = box->upper[d];
        }
    }
    return b;
}

void insertTree(Particle *p, TreeNode *t){
    Box sonbox;
    int b=sonNumber(&t->box, &sonbox, p);

    if(t->son[b] == NULL){
        if(t->node == particle){
            Particle p2 = t->p;
            t->son[b] = (TreeNode*)calloc(1,sizeof(TreeNode));
            t->son[b]->box = sonbox;
            t->son[b]->p = *p;

            t->son[b]->node = particle;
            t->node = pseudoParticle;

            insertTree(&p2, t);
        }else{
            t->son[b] = (TreeNode*)calloc(1,sizeof(TreeNode));
            t->son[b]->box = sonbox;
            t->son[b]->p = *p;

            t->son[b]->node = particle;
        }
    }else{
        insertTree(p, t->son[b]);
    }
}

