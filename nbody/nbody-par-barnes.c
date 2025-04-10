/*  
    Parallel N-Body simulation code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <mpi.h>


#define GRAVITY     1.1
#define FRICTION    0.01
#define MAXBODIES   10000
#define DELTA_T     (0.025/5000)
#define BOUNCE      -0.9
#define SEED        27102015


struct bodyType {
    double x[2];        /* Old and new X-axis coordinates */
    double y[2];        /* Old and new Y-axis coordinates */
    double xf;          /* force along X-axis */
    double yf;          /* force along Y-axis */
    double xv;          /* velocity along X-axis */
    double yv;          /* velocity along Y-axis */
    double mass;        /* Mass of the body */
    double radius;      /* width (derived from mass) */
};


struct world {
    struct bodyType bodies[MAXBODIES];
    int                 bodyCt;
    int                 old;    // Flips between 0 and 1

    /*  Dimensions of space (very finite, ain't it?) */
    int                 xdim;
    int                 ydim;
};


/* Quad Tree, 1 node has 4 child*/
struct quadTree {

    /* Child node, represent space partition*/
    /*
        [ TL ][ TR ]
        [ BL ][ BR ]
    */
    struct quadTree * child_TL;
    struct quadTree * child_TR;
    struct quadTree * child_BL;
    struct quadTree * child_BR;

    /* Tree node content */
    struct bodyType * body;

    /* Related Space Dimensions*/
    double xdim;
    double ydim;
    double xCenter;
    double yCenter;

    /* If node already has body, need to be sub divided */
    int has_body;

    /* If current tree is sub divied*/
    int divided;
    
    /* If node is parent node*/
    int parent;
};

/* Linked List to link all leave node */
struct linkedList
{
  struct bodyType * body;
  struct linkedList* next;

};

/*  Macros to hide memory layout
*/
#define X(w, B)        (w)->bodies[B].x[(w)->old]
#define XN(w, B)       (w)->bodies[B].x[(w)->old^1]
#define Y(w, B)        (w)->bodies[B].y[(w)->old]
#define YN(w, B)       (w)->bodies[B].y[(w)->old^1]
#define XF(w, B)       (w)->bodies[B].xf
#define YF(w, B)       (w)->bodies[B].yf
#define XV(w, B)       (w)->bodies[B].xv
#define YV(w, B)       (w)->bodies[B].yv
#define R(w, B)        (w)->bodies[B].radius
#define M(w, B)        (w)->bodies[B].mass

/* Visit linked list*/
static void
visit_list(struct linkedList ** head,int old){
    struct linkedList * cursor = * head;
    printf("list:\n");
    while(cursor  != NULL){
        printf("X: %f, Y: %f, Mass: %f\n",cursor->body->x[old],cursor->body->y[old],cursor->body->mass);
        cursor  = cursor -> next;
    }
}

/* Linked list insert*/
static void
insert_list(struct linkedList ** head, struct bodyType * body){
    struct linkedList * listNode =(struct linkedList *) malloc(sizeof(struct linkedList));
    struct linkedList * cursor = * head;
    listNode -> body = body;
    listNode -> next = NULL;
    if(*head == NULL){
        *head = listNode;
        return;
    }
    while(cursor->next != NULL)cursor = cursor -> next;
    cursor->next = listNode;
}

static void
init_node(struct quadTree * node,double xl,double xu,double yl,double yu){
    node->xdim = xu-xl;
    node->ydim = yu-yl;
    node->xCenter = (xl+xu) / 2;
    node->yCenter = (yl+yu) / 2;
    node->has_body = 0;
    node->divided = 0;
    node->parent = 0;
    node->child_BL = NULL;
    node->child_BR = NULL;
    node->child_TL = NULL;
    node->child_TR = NULL;
}

static void
print_tree(struct quadTree * root,int old,int depth)
{

    printf("Depth : %d\n",depth);
    printf("Parent  X: %f, Y: %f, Mass:%f\n",root->body->x[old],root->body->y[old],root->body->mass);
    if(root->child_TL != NULL){
        printf("TL child  X: %f, Y: %f, Mass:%f\n",root->child_TL->body->x[old],root->child_TL->body->y[old],root->child_TL->body->mass);
    }
    if(root->child_TR != NULL){
        printf("TR child  X: %f, Y: %f, Mass:%f\n",root->child_TR->body->x[old],root->child_TR->body->y[old],root->child_TR->body->mass);
    }
    if(root->child_BR != NULL){
        printf("BR child  X: %f, Y: %f, Mass:%f\n",root->child_BR->body->x[old],root->child_BR->body->y[old],root->child_BR->body->mass);
    }
    if(root->child_BL != NULL){
        printf("BL child  X: %f, Y: %f, Mass:%f\n",root->child_BL->body->x[old],root->child_BL->body->y[old],root->child_BL->body->mass);
    }

    if(root->child_TL != NULL){
        print_tree(root->child_TL,old,depth+1);
    }
    if(root->child_TR != NULL){
        print_tree(root->child_TR,old,depth+1);
    }
    if(root->child_BR != NULL){
        print_tree(root->child_BR,old,depth+1);
    }
    if(root->child_BL != NULL){
        print_tree(root->child_BL,old,depth+1);
    }
}

/* Insert new quad tree node */
static void
insert_tree(struct quadTree * root, struct bodyType * body,int old)
{

    if(root->has_body){
        if(!root->divided){

            /* 
                New coming body and old body stored in node,
                will be assign to correspond space partitioned 
                child node.
            */

            /* For new coming body */
            struct quadTree * newNode = malloc(sizeof(struct quadTree));
            struct  bodyType * newBody = malloc(sizeof(struct bodyType));
            memcpy(newBody,body,sizeof(struct bodyType));
            
            /* Space partition, insert new coming body */
            /* Body on boundry will be allocate to clock wise partition */
            /* Body on center will be allocate to TOP LEFT partition*/
            if( body->x[old] < root ->xCenter && ( body->y[old] > root ->yCenter || body->y[old] == root ->yCenter) ){
                root->child_TL = newNode;
                init_node(root->child_TL,   root ->xCenter - (root -> xdim)/2,
                                            root ->xCenter,
                                            root ->yCenter,
                                            root ->yCenter + (root -> ydim)/2);
                insert_tree(root->child_TL,newBody,old);
            }
            else if( ( body->x[old] > root ->xCenter || body->x[old] == root ->xCenter) && body->y[old] > root ->yCenter){
                root->child_TR = newNode;
                init_node(root->child_TR,   root ->xCenter,
                                            root ->xCenter + (root -> xdim)/2,
                                            root ->yCenter,
                                            root ->yCenter + (root -> ydim)/2);
                insert_tree(root->child_TR,newBody,old);
            }
            else if( ( body->x[old] < root ->xCenter || body->x[old] == root ->xCenter) && body->y[old] < root ->yCenter){
                root->child_BL = newNode;
                init_node(root->child_BL,   root ->xCenter - (root -> xdim)/2,
                                            root ->xCenter,
                                            root ->yCenter - (root -> ydim)/2,
                                            root ->yCenter);
                insert_tree(root->child_BL,newBody,old);
            }
            else if(body->x[old] > root ->xCenter && ( body->y[old] < root ->yCenter ||  body->y[old] == root ->yCenter) ){
                root->child_BR = newNode;
                init_node(root->child_BR,   root ->xCenter,
                                            root ->xCenter + (root -> xdim)/2,
                                            root ->yCenter - (root -> ydim)/2,
                                            root ->yCenter);
                insert_tree(root->child_BR,newBody,old);
            }
            else if(body->x[old] == root ->xCenter && body->y[old] == root ->yCenter ){
                root->child_TL = newNode;
                init_node(root->child_TL,   root ->xCenter - (root -> xdim)/2,
                                            root ->xCenter,
                                            root ->yCenter,
                                            root ->yCenter + (root -> ydim)/2);
                insert_tree(root->child_TL,newBody,old);
            }
            
            
            struct  bodyType * oldBody = malloc(sizeof(struct bodyType));
            memcpy(oldBody,root->body,sizeof(struct bodyType));

            /* Space partition, insert stored old body to child leave */
            if( (root->body)->x[old] < root ->xCenter && ( ((root->body)->y[old] > root ->yCenter) || ((root->body)->y[old] == root ->yCenter))){
                if(root->child_TL != NULL){
                    insert_tree(root->child_TL,oldBody,old);
                }else{
                    struct quadTree * oldNode = malloc(sizeof(struct quadTree));
                    root->child_TL = oldNode;
                    init_node(root->child_TL,   root ->xCenter - (root -> xdim)/2,
                                                root ->xCenter,
                                                root ->yCenter,
                                                root ->yCenter + (root -> ydim)/2);
                    insert_tree(root->child_TL,oldBody,old);
                }
            }
            else if( ( (root->body)->x[old] > root ->xCenter || (root->body)->x[old] == root ->xCenter)&& (root->body)->y[old] > root ->yCenter){
                if(root->child_TR != NULL){
                    insert_tree(root->child_TR,oldBody,old);
                }else{
                    struct quadTree * oldNode = malloc(sizeof(struct quadTree));
                    root->child_TR = oldNode;
                    init_node(root->child_TR,   root ->xCenter,
                                                root ->xCenter + (root -> xdim)/2,
                                                root ->yCenter,
                                                root ->yCenter + (root -> ydim)/2);
                    insert_tree(root->child_TR,oldBody,old);
                }
            }
            else if( ( (root->body)->x[old] < root ->xCenter || (root->body)->x[old] == root ->xCenter)&& (root->body)->y[old] < root ->yCenter){
                if(root->child_BL != NULL){
                    insert_tree(root->child_BL,oldBody,old);
                }else{
                    struct quadTree * oldNode = malloc(sizeof(struct quadTree));
                    root->child_BL = oldNode;
                    init_node(root->child_BL,   root ->xCenter - (root -> xdim)/2,
                                                root ->xCenter,
                                                 root ->yCenter - (root -> ydim)/2,
                                                root ->yCenter);
                    insert_tree(root->child_BL,oldBody,old);
                }    
            }
            else if((root->body)->x[old] > root ->xCenter && ( (root->body)->y[old] < root ->yCenter || (root->body)->y[old] == root ->yCenter ) ){
                if(root->child_BR != NULL){
                    insert_tree(root->child_BR,oldBody,old);
                }else{
                    struct quadTree * oldNode = malloc(sizeof(struct quadTree));
                    root->child_BR = oldNode;
                    init_node(root->child_BR,   root ->xCenter,
                                                root ->xCenter + (root -> xdim)/2,
                                                root ->yCenter - (root -> ydim)/2,
                                                root ->yCenter);
                    insert_tree(root->child_BR,oldBody,old);
                }
            }
            else if((root->body)->x[old] == root ->xCenter && (root->body)->y[old] == root ->yCenter ){
                if(root->child_TL != NULL){
                    insert_tree(root->child_TL,oldBody,old);
                }else{
                    struct quadTree * oldNode = malloc(sizeof(struct quadTree));
                    root->child_TL = oldNode;
                    init_node(root->child_TL,   root ->xCenter - (root -> xdim)/2,
                                                root ->xCenter,
                                                root ->yCenter,
                                                root ->yCenter + (root -> ydim)/2);
                    insert_tree(root->child_TL,oldBody,old);
                }
            }

            root->divided = 1;
            root->parent = 1;

            /*
            struct  bodyType * parentBody = malloc(sizeof(struct bodyType));
            memset(parentBody, 0, sizeof(struct bodyType));
            root->body = parentBody;
            */

        }else{

            /*  
                If the tree node is already divied,
                insert new coming node into child node,
                check if child node is empty before insertion
            */

            if( body->x[old] < root ->xCenter && ( body->y[old] > root ->yCenter || body->y[old] == root ->yCenter) ){
                if(root->child_TL != NULL){
                    insert_tree(root->child_TL,body,old);
                }else{
                    struct quadTree * newNode = malloc(sizeof(struct quadTree));
                    struct  bodyType * newBody = malloc(sizeof(struct bodyType));
                    memcpy(newBody,body,sizeof(struct bodyType));
                    root->child_TL = newNode;
                    init_node(root->child_TL,   root ->xCenter - (root -> xdim)/2,
                                                root ->xCenter,
                                                root ->yCenter,
                                                root ->yCenter + (root -> ydim)/2);
                    insert_tree(root->child_TL,newBody,old);
                }
            }
            else if( ( body->x[old] > root ->xCenter || body->x[old] == root ->xCenter) && body->y[old] > root ->yCenter){
                if(root->child_TR != NULL){
                    insert_tree(root->child_TR,body,old);
                }else{
                    struct quadTree * newNode = malloc(sizeof(struct quadTree));
                    struct  bodyType * newBody = malloc(sizeof(struct bodyType));
                    memcpy(newBody,body,sizeof(struct bodyType));
                    root->child_TR = newNode;
                    init_node(root->child_TR,   root ->xCenter,
                                                root ->xCenter + (root -> xdim)/2,
                                                root ->yCenter,
                                                root ->yCenter + (root -> ydim)/2);
                    insert_tree(root->child_TR,newBody,old);
                }

            }
            else if( ( body->x[old] < root ->xCenter || body->x[old] == root ->xCenter) && body->y[old] < root ->yCenter){
                if(root->child_BL != NULL){
                    insert_tree(root->child_BL,body,old);
                }else{
                    struct quadTree * newNode = malloc(sizeof(struct quadTree));
                    struct  bodyType * newBody = malloc(sizeof(struct bodyType));
                    memcpy(newBody,body,sizeof(struct bodyType));
                    root->child_BL = newNode;
                    init_node(root->child_BL,   root ->xCenter - (root -> xdim)/2,
                                                root ->xCenter,
                                                root ->yCenter - (root -> ydim)/2,
                                                root ->yCenter);
                    insert_tree(root->child_BL,newBody,old);
                }
                
            }
            else if(body->x[old] > root ->xCenter && ( body->y[old] < root ->yCenter ||  body->y[old] == root ->yCenter) ){
                if(root->child_BR != NULL){
                    insert_tree(root->child_BR,body,old);
                }else{
                    struct quadTree * newNode = malloc(sizeof(struct quadTree));
                    struct  bodyType * newBody = malloc(sizeof(struct bodyType));
                    memcpy(newBody,body,sizeof(struct bodyType));
                    root->child_BR = newNode;
                    init_node(root->child_BR,   root ->xCenter,
                                                root ->xCenter + (root -> xdim)/2,
                                                root ->yCenter - (root -> ydim)/2,
                                                root ->yCenter);
                    insert_tree(root->child_BR,newBody,old);
                }
            }
            else if(body->x[old] == root ->xCenter && body->y[old] == root ->yCenter ){
                if(root->child_TL != NULL){
                    insert_tree(root->child_TL,body,old);
                }else{
                    struct quadTree * newNode = malloc(sizeof(struct quadTree));
                    struct  bodyType * newBody = malloc(sizeof(struct bodyType));
                    memcpy(newBody,body,sizeof(struct bodyType));
                    root->child_TL = newNode;
                    init_node(root->child_TL,   root ->xCenter - (root -> xdim)/2,
                                                root ->xCenter,
                                                root ->yCenter,
                                                root ->yCenter + (root -> ydim)/2);
                    insert_tree(root->child_TL,newBody,old);
                }
            }
            
        }
    }else{
        root->body = body;
        root->has_body = 1;
    }
    
}

/* Build Quadtree*/
static void
build_tree(struct world *world,struct quadTree * root)
{
    init_node(root,0,world->xdim,0,world->ydim);

    for(int b = 0; b < world->bodyCt; ++b){

        //printf("x: %f, y: %f \n",world->bodies[b].x[world->old],world->bodies[b].y[world->old]);
        insert_tree(root,&(world->bodies[b]),world->old);
    }
}

/* Post order traversal to compute center of mass */
/* Visting child starting from TL, clock-wise*/
/* Save visited node on linked list */
static void
compute_center_mass(struct quadTree * root,int old, struct linkedList ** head){
    if(root == NULL) return;
    else{
        if(root->parent){
            compute_center_mass(root->child_TL,old,head);
            compute_center_mass(root->child_TR,old,head);
            compute_center_mass(root->child_BR,old,head);
            compute_center_mass(root->child_BL,old,head);


            int mass_sum = 0;
            int x_sum = 0;
            int y_sum = 0;

            if(root->child_TL != NULL){
                mass_sum += root->child_TL->body->mass;
                x_sum += root->child_TL->body->x[old] * root->child_TL->body->mass;
                y_sum += root->child_TL->body->y[old] * root->child_TL->body->mass;
                //printf("TL child  X: %f, Y: %f,     Mass:%f \n",root->child_TL->body->x[old],root->child_TL->body->y[old],root->child_TL->body->mass);
            }
            if(root->child_TR != NULL){
                mass_sum += root->child_TR->body->mass;
                x_sum += root->child_TR->body->x[old] * root->child_TR->body->mass;
                y_sum += root->child_TR->body->y[old] * root->child_TR->body->mass;
                //printf("TR child  X: %f, Y: %f,     Mass:%f \n",root->child_TR->body->x[old],root->child_TR->body->y[old],root->child_TR->body->mass);
            }
            if(root->child_BR != NULL){
                mass_sum += root->child_BR->body->mass;
                x_sum += root->child_BR->body->x[old] * root->child_BR->body->mass;
                y_sum += root->child_BR->body->y[old] * root->child_BR->body->mass;
                //printf("BR child  X: %f, Y: %f,     Mass:%f \n",root->child_BR->body->x[old],root->child_BR->body->y[old],root->child_BR->body->mass);
            }
            if(root->child_BL != NULL){
                mass_sum += root->child_BL->body->mass;
                x_sum += root->child_BL->body->x[old] * root->child_BL->body->mass;
                y_sum += root->child_BL->body->y[old] * root->child_BL->body->mass;
                //printf("BL child  X: %f, Y: %f,     Mass:%f \n",root->child_BL->body->x[old],root->child_BL->body->y[old],root->child_BL->body->mass);
            }
            root->body->mass = mass_sum;
            root->body->x[old] = x_sum/mass_sum;
            root->body->y[old] = y_sum/mass_sum;
            

        }else{
            insert_list(head,root->body);
        }
    }

}


/* Clear accumlated forces on each astro body*/
static void
clear_forces(struct world *world)
{
    int b;
    for (b = 0; b < world->bodyCt; ++b) {
        YF(world, b) = XF(world, b) = 0;
    }
}

static void 
compute_forces(struct world *world)
{
    int b, c;

    /* Incrementally accumulate forces from each body pair,
       skipping force of body on itself (c == b)
    */
    for (b = 0; b < world->bodyCt; ++b) {
        for (c = b + 1; c < world->bodyCt; ++c) {
            double dx = X(world, c) - X(world, b);
            double dy = Y(world, c) - Y(world, b);
            double angle = atan2(dy, dx);
            double dsqr = dx*dx + dy*dy;
            double mindist = R(world, b) + R(world, c);
            double mindsqr = mindist*mindist;
            // If distance between x and y is bigger than their radius
            double forced = ((dsqr < mindsqr) ? mindsqr : dsqr);
            double force = M(world, b) * M(world, c) * GRAVITY / forced;
            double xf = force * cos(angle);
            double yf = force * sin(angle);

            /* Slightly sneaky...
               force of b on c is negative of c on b;
            */
            XF(world, b) += xf;
            YF(world, b) += yf;
            XF(world, c) -= xf;
            YF(world, c) -= yf;
        }
    }
}

/* Compute velocity from accumlated foreces */
static void
compute_velocities(struct world *world)
{
    int b;

    for (b = 0; b < world->bodyCt; ++b) {
        double xv = XV(world, b);
        double yv = YV(world, b);
        double force = sqrt(xv*xv + yv*yv) * FRICTION;
        double angle = atan2(yv, xv);
        double xf = XF(world, b) - (force * cos(angle));
        double yf = YF(world, b) - (force * sin(angle));

        XV(world, b) += (xf / M(world, b)) * DELTA_T;
        YV(world, b) += (yf / M(world, b)) * DELTA_T;
    }
}

/* Compute astro body position from velocity and time interval */
static void
compute_positions(struct world *world)
{
    int b;

    for (b = 0; b < world->bodyCt; ++b) {
        double xn = X(world, b) + XV(world, b) * DELTA_T;
        double yn = Y(world, b) + YV(world, b) * DELTA_T;

        /* Bounce off image "walls" */
        if (xn < 0) {
            xn = 0;
            XV(world, b) = -XV(world, b);
        } else if (xn >= world->xdim) {
            xn = world->xdim - 1;
            XV(world, b) = -XV(world, b);
        }
        if (yn < 0) {
            yn = 0;
            YV(world, b) = -YV(world, b);
        } else if (yn >= world->ydim) {
            yn = world->ydim - 1;
            YV(world, b) = -YV(world, b);
        }

        /* Update position */
        XN(world, b) = xn;
        YN(world, b) = yn;
    }
}

/*  */

/*  Graphic output stuff...
 */

#include <fcntl.h>
#include <sys/mman.h>

struct filemap {
    int            fd;
    off_t          fsize;
    void          *map;
    unsigned char *image;
};


static void
filemap_close(struct filemap *filemap)
{
    if (filemap->fd == -1) {
        return;
    }
    close(filemap->fd);
    if (filemap->map == MAP_FAILED) {
        return;
    }
    munmap(filemap->map, filemap->fsize);
}

static unsigned char *
Eat_Space(unsigned char *p)
{
    while ((*p == ' ') ||
           (*p == '\t') ||
           (*p == '\n') ||
           (*p == '\r') ||
           (*p == '#')) {
        if (*p == '#') {
            while (*(++p) != '\n') {
                // skip until EOL
            }
        }
        ++p;
    }

    return p;
}

static unsigned char *
Get_Number(unsigned char *p, int *n)
{
    p = Eat_Space(p);  /* Eat white space and comments */

    int charval = *p;
    if ((charval < '0') || (charval > '9')) {
        errno = EPROTO;
        return NULL;
    }

    *n = (charval - '0');
    charval = *(++p);
    while ((charval >= '0') && (charval <= '9')) {
        *n *= 10;
        *n += (charval - '0');
        charval = *(++p);
    }

    return p;
}

static int
map_P6(const char *filename, int *xdim, int *ydim, struct filemap *filemap)
{
    /* The following is a fast and sloppy way to
       map a color raw PPM (P6) image file
    */
    int maxval;
    unsigned char *p;

    /* First, open the file... */
    if ((filemap->fd = open(filename, O_RDWR)) < 0) {
        goto ppm_abort;
    }

    /* Read size and map the whole file... */
    filemap->fsize = lseek(filemap->fd, (off_t)0, SEEK_END);
    filemap->map = mmap(0,                      // Put it anywhere
                        filemap->fsize,         // Map the whole file
                        (PROT_READ|PROT_WRITE), // Read/write
                        MAP_SHARED,             // Not just for me
                        filemap->fd,            // The file
                        0);                     // Right from the start
    if (filemap->map == MAP_FAILED) {
        goto ppm_abort;
    }

    /* File should now be mapped; read magic value */
    p = filemap->map;
    if (*(p++) != 'P') goto ppm_abort;
    switch (*(p++)) {
    case '6':
        break;
    default:
        errno = EPROTO;
        goto ppm_abort;
    }

    p = Get_Number(p, xdim);            // Get image width */
    if (p == NULL) goto ppm_abort;
    p = Get_Number(p, ydim);            // Get image width */
    if (p == NULL) goto ppm_abort;
    p = Get_Number(p, &maxval);         // Get image max value */
    if (p == NULL) goto ppm_abort;

    /* Should be 8-bit binary after one whitespace char... */
    if (maxval > 255) {
        goto ppm_abort;
    }
    if ((*p != ' ') &&
        (*p != '\t') &&
        (*p != '\n') &&
        (*p != '\r')) {
        errno = EPROTO;
        goto ppm_abort;
    }

    /* Here we are... next byte begins the 24-bit data */
    filemap->image = p + 1;

    return 0;

ppm_abort:
    filemap_close(filemap);

    return -1;
}

static inline void
color(const struct world *world, unsigned char *image, int x, int y, int b)
{
    unsigned char *p = image + (3 * (x + (y * world->xdim)));
    int tint = ((0xfff * (b + 1)) / (world->bodyCt + 2));

    p[0] = (tint & 0xf) << 4;
    p[1] = (tint & 0xf0);
    p[2] = (tint & 0xf00) >> 4;
}

static inline void
black(const struct world *world, unsigned char *image, int x, int y)
{
    unsigned char *p = image + (3 * (x + (y * world->xdim)));

    p[2] = (p[1] = (p[0] = 0));
}

static void
display(const struct world *world, unsigned char *image)
{
    double i, j;
    int b;

    /* For each pixel */
    for (j = 0; j < world->ydim; ++j) {
        for (i = 0; i < world->xdim; ++i) {
            /* Find the first body covering here */
            for (b = 0; b < world->bodyCt; ++b) {
                double dy = Y(world, b) - j;
                double dx = X(world, b) - i;
                double d = sqrt(dx*dx + dy*dy);

                if (d <= R(world, b)+0.5) {
                    /* This is it */
                    color(world, image, i, j, b);
                    goto colored;
                }
            }

            /* No object -- empty space */
            black(world, image, i, j);

colored:        ;
        }
    }
}

static void
print(struct world *world)
{
    int b;

    for (b = 0; b < world->bodyCt; ++b) {
        printf("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
               X(world, b), Y(world, b), XF(world, b), YF(world, b), XV(world, b), YV(world, b));
    }
}

static long long
nr_flops(int n, int steps) {
  long long nr_flops = 0;
  // compute forces
  nr_flops += 20 * (n * (n-1) / 2);
  // compute velocities
  nr_flops += 18 * n;
  // compute positions
  nr_flops += 4 * n;

  nr_flops *= steps;

  return nr_flops;
}
    

/*  Main program...
*/

int
main(int argc, char **argv)
{
    unsigned int lastup = 0;
    unsigned int secsup;
    int b;
    int steps;
    double rtime;
    struct timeval start;
    struct timeval end;
    struct filemap image_map;

    // MPI Initlization 
    MPI_Init(NULL, NULL);
    int MPI_rank;
    int MPI_world;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);   // GET current node rank
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_world);  // GET total node size


    struct world *world = calloc(1, sizeof *world);
    if (world == NULL) {
        fprintf(stderr, "Cannot calloc(world)\n");
        exit(1);
    }

    /* Get Parameters */
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_bodies secs_per_update ppm_output_file steps\n",
                argv[0]);
        exit(1);
    }
    if ((world->bodyCt = atol(argv[1])) > MAXBODIES ) {
        fprintf(stderr, "Using only %d bodies...\n", MAXBODIES);
        world->bodyCt = MAXBODIES;
    } else if (world->bodyCt < 2) {
        fprintf(stderr, "Using two bodies...\n");
        world->bodyCt = 2;
    }
    secsup = atoi(argv[2]);
    if (map_P6(argv[3], &world->xdim, &world->ydim, &image_map) == -1) {
        fprintf(stderr, "Cannot read %s: %s\n", argv[3], strerror(errno));
        exit(1);
    }
    steps = atoi(argv[4]);

    fprintf(stderr, "Running N-body with %i bodies and %i steps\n", world->bodyCt, steps);

    /* Initialize simulation data */
    srand(SEED);
    for (b = 0; b < world->bodyCt; ++b) {
        X(world, b) = (rand() % world->xdim);
        Y(world, b) = (rand() % world->ydim);
        R(world, b) = 1 + ((b*b + 1.0) * sqrt(1.0 * ((world->xdim * world->xdim) + (world->ydim * world->ydim)))) /
                (25.0 * (world->bodyCt*world->bodyCt + 1.0));
        M(world, b) = R(world, b) * R(world, b) * R(world, b);
        XV(world, b) = ((rand() % 20000) - 10000) / 2000.0;
        YV(world, b) = ((rand() % 20000) - 10000) / 2000.0;
    }

    if (gettimeofday(&start, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    /* Build Tree */
    if(MPI_rank == 0){
        struct quadTree * treeRoot = malloc(sizeof(struct quadTree));
        memset(treeRoot, 0, sizeof(struct quadTree)); 
        build_tree(world,treeRoot);
        //print_tree(treeRoot,world->old,0);
        struct linkedList * listHead = NULL;
        compute_center_mass(treeRoot,world->old,&listHead);
        visit_list(&listHead,world->old);
    }

    /* Main Loop */
    for (int step = 0; step < steps; step++) {
        clear_forces(world);
        compute_forces(world);
        compute_velocities(world);
        compute_positions(world);

        /* Flip old & new coordinates */
        world->old ^= 1;

        /*Time for a display update?*/ 
        if (secsup > 0 && (time(0) - lastup) > secsup) {
            display(world, image_map.image);
            msync(image_map.map, image_map.fsize, MS_SYNC); /* Force write */
            lastup = time(0);
        }
    }

    if (gettimeofday(&end, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    rtime = (end.tv_sec + (end.tv_usec / 1000000.0)) - 
                (start.tv_sec + (start.tv_usec / 1000000.0));


    print(world);

    fprintf(stderr, "\nN-body took: %.3f seconds\n", rtime);
    fprintf(stderr, "Performance N-body: %.2f GFLOPS\n", nr_flops(world->bodyCt, steps) / 1e9 / rtime);

    filemap_close(&image_map);

    free(world);

    MPI_Finalize();

    return 0;
}
