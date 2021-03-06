struct devent {
	double time;
	int popi;
	int popj;
	double paramv;
	double **mat ;
	char detype ;
	struct devent *nextde;
	} ;
struct c_params {
	int npop;
	int nsam;
        int **config ;
	double **mig_mat;
	double r;
	int nsites;
	double f;
	double track_len;
	double *size;
	double *alphag;
        int selspot ;
		int nrepsfreq;
        int nfreqder ;
        int maxpop ;
        double **infreqder ;
        double *tfreqder ;
	struct devent *deventlist ;
	} ;
struct m_params {
	 double theta;
	int segsitesin;
	int treeflag;
	int timeflag;
	int mfreq;
	 } ;
struct params { 
	struct c_params cp;
	struct m_params mp;
        FILE *selfile ;
        int commandlineseedflag ;
		int selsitepolyflag ;
		int output_precision;
	};
	
struct node{
	int abv;
	int ndes;
	float time;
};


void ordran(int n, double pbuf[]);
void ranvec(int n, double pbuf[]);
void order(int n, double pbuf[]);

void biggerlist(int nsam,  char **list );
int poisso(double u);
void locate(int n,double beg, double len,double *ptr);
void mnmial(int n, int nclass, double p[], int rv[]);
void usage();
int tdesn(struct node *ptree, int tip, int node );
int pick2(int n, int *i, int *j);
int xover(int nsam,int ic, int is);
int links(int c);


