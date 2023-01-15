/*  
    N-Body simulation code.
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

struct job
{
    int bodyPar[2];
};

struct jobList
{

    struct job * job;
    struct jobList * next;
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

/* Append linked list*/
static struct jobList **
append_list(struct jobList ** jobList,int x1,int x2)
{
    struct jobList * jobNode = (struct jobList *)malloc(sizeof(struct jobList));
    struct job * job= (struct job *)malloc(sizeof(struct job));
    job->bodyPar[0] = x1;
    job->bodyPar[1] = x2;
    jobNode->job = job;
    jobNode->next = NULL;
    *jobList = jobNode;
    return &(jobNode->next);
}

/* Generate job for different node*/
static int
generate_job(struct world *world,struct jobList ** jobList, int MPI_world,int MPI_rank)
{
    int amount = 0;
    /* Jobs are seprated into sections */
    int sec_begin,sec_end;
    int jobIndex = 0;

    struct jobList ** cursor = jobList;
    for(int i = 0; i < world->bodyCt; i++)
    {
        for (int j = i+1; j < world->bodyCt; j++)
        {
            amount += 1;            
        }
    }
    sec_begin = MPI_rank * (int)(amount / MPI_world);
    sec_end = sec_begin + amount / MPI_world;
    if(MPI_rank == MPI_world - 1)sec_end = amount;
    
    for(int i = 0; i < world->bodyCt; i++)
    {
        for (int j = i+1; j < world->bodyCt; j++)
        {           
            /* If job index is within section, add to job list*/
            if( (jobIndex > sec_begin || jobIndex == sec_begin)&&(jobIndex < sec_end || jobIndex == sec_end) ){
                cursor = append_list(cursor,i,j);
            }
            jobIndex += 1;
            if (jobIndex > sec_end)
            {
                break;
            }
        }
        if (jobIndex > sec_end)
        {
            break;
        }
    }
    
    return amount;
}

static void
clear_forces(struct world *world)
{
    int b;

    /* Clear force accumulation variables */
    for (b = 0; b < world->bodyCt; ++b) {
        YF(world, b) = XF(world, b) = 0;
    }
}

static void
compute_forces(struct world *world,struct jobList ** jobList)
{
    int b, c;

    //list_traversal(jobList);
    /* Incrementally accumulate forces from each body pair,
       skipping force of body on itself (c == b)
    */
   struct jobList * cursor = * jobList;
    while(cursor  != NULL){
        b = cursor->job->bodyPar[0];
        c = cursor->job->bodyPar[1];
        double dx = X(world, c) - X(world, b);
        double dy = Y(world, c) - Y(world, b);
        double angle = atan2(dy, dx);
        double dsqr = dx*dx + dy*dy;
        double mindist = R(world, b) + R(world, c);
        double mindsqr = mindist*mindist;
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
        cursor  = cursor -> next;
    }
    
}


/* Version 1.0 Gather*/
/* Possible : sequential gather to hide latency*/
static void
merge_force(double * world_X_force,double * world_Y_force,
            double * parallel_X_force,double * parallel_Y_force,
            int MPI_world,struct world *world,int rank)
{
    
    for(int i =0;i<world->bodyCt;i++){
        world_X_force[i] = world->bodies[i].xf;
        world_Y_force[i] = world->bodies[i].yf;
    }
    
    MPI_Gather(world_X_force,world->bodyCt,MPI_DOUBLE,parallel_X_force,world->bodyCt,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(world_Y_force,world->bodyCt,MPI_DOUBLE,parallel_Y_force,world->bodyCt,MPI_DOUBLE,0,MPI_COMM_WORLD);

    int off_set = 0;
    if(rank == 0){
        /* Merge force from different world to root world */
        for(int i = 1;i < MPI_world;i++){
            off_set = i * world->bodyCt;
            for(int j = 0; j < world->bodyCt; j++){
                world_X_force[j] +=  parallel_X_force[j + off_set];
                world_Y_force[j] +=  parallel_Y_force[j + off_set];
            }
        }
        /* Update root world bodies  */
        for(int j = 0; j < world->bodyCt; j++){
                XF(world, j) = world_X_force[j];
                YF(world, j) = world_Y_force[j];
        }
    }
    
}


/* Version 1.0 Broadcast*/
/* Root node update other nodes*/
static void
update_world(   double * world_buf_X_1,double * world_buf_Y_1,
                double * world_buf_X_2,double * world_buf_Y_2,
                struct world *world,int rank)
{
    /* root store computed data into send buffer */
    if(rank == 0){
        for(int i =0;i<world->bodyCt;i++){
            world_buf_X_1[i] = XV(world, i);
            world_buf_Y_1[i] = YV(world, i);
            world_buf_X_2[i] = XN(world, i);
            world_buf_Y_2[i] = YN(world, i);
        }
    }   
    
    MPI_Bcast(world_buf_X_1,world->bodyCt,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(world_buf_Y_1,world->bodyCt,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(world_buf_X_2,world->bodyCt,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(world_buf_Y_2,world->bodyCt,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for(int i =0;i<world->bodyCt;i++){
            XV(world, i) = world_buf_X_1[i];
            YV(world, i) = world_buf_Y_1[i];
            XN(world, i) = world_buf_X_2[i];
            YN(world, i) = world_buf_Y_2[i];
        }
    
}

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

    struct world * world = calloc(1, sizeof *world);

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

    if(MPI_rank == 0){

    print(world);
    printf("***************************\n");
    }

    /* Store body XY forces for parallel communication */
    double * world_X_force = malloc(world->bodyCt*(sizeof(double)));
    double * world_Y_force = malloc(world->bodyCt*(sizeof(double)));


    double * world_buf_X_2 = malloc(world->bodyCt*(sizeof(double)));
    double * world_buf_Y_2 = malloc(world->bodyCt*(sizeof(double)));

    double * parallel_X_force = NULL;
    double * parallel_Y_force = NULL;

    /* If root node, create space to store forces from different node*/
    if(MPI_rank == 0){
        parallel_X_force = malloc((MPI_world * world->bodyCt) * (sizeof(double)));
        parallel_Y_force = malloc((MPI_world * world->bodyCt) * (sizeof(double)));
        if (parallel_X_force == NULL || parallel_Y_force == NULL) {
        fprintf(stderr, "Cannot calloc(parallel foreces)\n");
        exit(1);
        }
    }


    if (gettimeofday(&start, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    
    struct jobList * jobList = NULL;
    generate_job(world,&jobList,MPI_world,MPI_rank);

    int count =steps;
    /* Main Loop */
    for (int step = 0; step < steps; step++) {


        clear_forces(world);
        compute_forces(world,&jobList);
        /* Root gather force data from other node */
        /* Merege force on bodies from different world */
        merge_force(world_X_force,world_Y_force,parallel_X_force,parallel_Y_force,MPI_world,world,MPI_rank);
        compute_velocities(world);

        /* Broadcast new computed velocities to other node */
        /*
        MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
               MPI_Comm comm)
        if(MPI_rank == 0 && count == 1){
            for(int i =0 ; i < world->bodyCt;i++)printf("X: %f, Y: %f\n",XV(world, i),YV(world, i));
            printf("********************************************\n");
            count = 0;
        }*/
        compute_positions(world);
        update_world(world_X_force,world_Y_force,world_buf_X_2,world_buf_Y_2,world,MPI_rank);
        count--;
       
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


    if(MPI_rank == 0){
        print(world);
        fprintf(stderr, "\nCore: %d", MPI_rank);
        fprintf(stderr, "\nN-body took: %.3f seconds\n", rtime);
        fprintf(stderr, "Performance N-body: %.2f GFLOPS\n", nr_flops(world->bodyCt, steps) / 1e9 / rtime);
    }

    filemap_close(&image_map);

    free(world);

    MPI_Finalize();

    return 0;
}
