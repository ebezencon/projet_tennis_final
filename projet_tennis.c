/*
Projet Tennis - Eliott Bezençon, SIE-BA3. Automne 2020.

Cette simulation à pour but d'évaluer la trajectoire et la faisabilité d'un coup dont
on connaît les conditions initiales. Les nombreux paramètres pris en compte ont pour 
objectif de simuler un coup aussi proche que possible de la réalité.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

//Initialisation des structures

struct vecteur {        //vecteur utilisé pour position,vitesse,acceleration,forces,impacts (coordonées des points d'impact)
    double x;
    double y;
    double z;
};

struct limites {        //Selon le type de coup que nous entreprenons, les limites du terrain seront différentes.
    char* type_de_coup;
    double x_max;
    double x_min;
    double y_max;
    double y_min;
    char* nom;          //ie carré de service ou terrain (nom des limites)
};

struct CST {        //structure contenant les constante propre au tennis
    double masse;
    double poids;
    double rho;      //densité air
    double r;        //rayon de la balle de tennis
    double s;        //coeff de magnus (pour lift&slice)
    double cf;       //coeff de frottement ds force trainée
    double wx;       //vitesse angulaire selon les trois axes, constantes durant tout l'envol
    double wy;
    double wz; 
    double e;        //coeff de restitution des vitesses (Energie cinétique) après collision avec le sol
};

struct settings {              //permettent de spécifier les aspects de la simulation sur lesquels on veut se concentrer.
    char* nom_fichier;
    char* frottement_air;
    char* calcule_faisabilite;
    int faisabilite_n;
    double dt;
    char* nom_fichier_impacts;
    char* type_de_coup;
    double rayon_de_similitude;
};

struct alea {                   //Partie aléatoire de la simulation qui nous donne des infos sur la faisabilité du coup.
    double ecart_type_vitesse;
    double ecart_type_angles;
    double ecart_type_w;
    double vent;
    double phi_vent;
    int dans_limites;           //on stocke le nombre de balles qui ont atteris dans les diff limites pour calculer les probas
    int dans_cercle;
    int dans_limites_et_dans_cercle;
};

//Fonctions

struct vecteur coordvitessecartesienne(double module_vitesse_initiale, double teta, double phi){
    //a partir d'un système de coordonée sphèrique on passe en coordonée cartésienne et renvoie le vecteur vitesse
    struct vecteur v;
    v.x=module_vitesse_initiale*sin(teta)*cos(phi);
    v.y=module_vitesse_initiale*sin(teta)*sin(phi);
    v.z=module_vitesse_initiale*cos(teta);
    return v;
}

void calcule_acceleration(struct vecteur * a, struct vecteur * forces, struct vecteur v, struct CST cst) {
    //Fonction qui calcul les accélerations à partir de la deuxième loi de newton
    double module_vitesse = sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
    //Calcul des forces selons les 3 axes
    forces->x = cst.cf*module_vitesse*v.x  +  cst.s*(cst.wy*v.z-cst.wz*v.y);              // force de trainee et force de magnus (produit vectoriel)
    forces->y = cst.cf*module_vitesse*v.y  +  cst.s*(cst.wz*v.x-cst.wx*v.z);              // force de trainee et force de magnus (produit vectoriel)
    forces->z = cst.cf*module_vitesse*v.z  +  cst.s*(cst.wx*v.y-cst.wy*v.x) + cst.poids;  // Poids, force trainee et force magnus
    //Deuxième loi newton : somme forces/masse = acceleration
    a->x=forces->x/cst.masse;
    a->y=forces->y/cst.masse;
    a->z=forces->z/cst.masse;
}

void calcule_position_et_vitesse(struct vecteur * x, struct vecteur * v, struct vecteur a,struct settings settings){
    //Développement limité de x(t+dt) d'ordre 2 pour la position
    x->x += v->x*settings.dt + 1/2*a.x*settings.dt*settings.dt;
    x->y += v->y*settings.dt + 1/2*a.y*settings.dt*settings.dt;
    x->z += v->z*settings.dt + 1/2*a.z*settings.dt*settings.dt;
    //Développement limité de v(t+dt) d'ordre 1 pour la vitesse
    v->x += a.x * settings.dt;
    v->y += a.y * settings.dt;
    v->z += a.z * settings.dt;
}

int filet(double intervale_filet,struct vecteur * x){   //intervalle filet = interval ni trop petit, ni trop grand dans lequel la balle passera une fois pendant son envol, x=vecteur position (avec trois coordonées)
    //Cette fonction retoune 1 si le coup part dans le filet, 0 sinon. (la hauteur du filet est 0.91m)
    //elle vérifie si la balle se situe d'abord dans l'interval très proche du filet et si c'est le cas, elle vérifie si l'on est au dessus de celui ci,à côté ou dedans.
    if (-intervale_filet/2 < x->y && x->y < intervale_filet/2 && x->z < 0.91 &&x->x>-5.49&&x->x<5.49) {
        return 1;
    } else {
        return 0;
    }
}   

struct limites determination_limites(char* type_de_coup){
    //cette fonction retourne une structure contenant les limites associée au type de coup rentré dans les settings, à l'intérieure desquelles la balle est IN dans les règles du tennis.
    struct limites limites;
    if (strcmp(type_de_coup,"service_gauche")==0){
        //rectangle de gauche
        limites.type_de_coup="service";
        limites.nom="carré de service";
        limites.x_max=0;
        limites.x_min=-4.12;
        limites.y_max=6.40;
        limites.y_min=0;
        return limites;
    } else if (strcmp(type_de_coup,"service_droite")==0){
        //rectangle de droite
        limites.type_de_coup="service";
        limites.nom="carré de service";
        limites.x_max=4.12;
        limites.x_min=0;
        limites.y_max=6.40;
        limites.y_min=0;
        return limites;
    } else {
        //limite du terrain en entier (sans les corridors)
        limites.type_de_coup="coup";
        limites.nom="terrain";
        limites.x_max=4.12;
        limites.x_min=-4.12;
        limites.y_max=11.89;
        limites.y_min=0;
        return limites;
    }
}

void simulation_sans_alea(struct vecteur * x, struct vecteur * v, struct vecteur a, struct vecteur forces,struct CST cst,struct settings* settings, struct vecteur* impact, struct limites limites){ 
    //Cette fonction est la simulation principale sur laquelle on relève tous les points de la trajectoire du coup pas impacté par les aléas
    //elle nous donne aussi des informations sur les propriétés du coups qu'on a fait ; si c'est un bon coup ou si il ratteri trop tot, si il a fini dans le filet ou en dehors des limites du terrain

    //On ouvre le fichier pour pouvoir y mettre les coordonnées de notre trajectoire
    FILE * fichier = fopen(settings->nom_fichier, "w");
    fprintf(fichier,"X,Y,Z\n");                                                  //première ligne avec X,Y,Z pour faciliter lecture fichier dans prog python
    
    double t = 0;
    double intervale_filet = v->y*settings->dt;                                  //interval dans lequel on vérifie si l'on passe le filet ou non (intervalle très petit) (on y sera de toute façon car vy baisse avec le temps.)
    double rebonds = 0;                                                          //on arrète la simulation après deux rebonds

    while (x->x <7 && x->x>-7 && x->y< 15 && x->y>-14 && rebonds<2) {            //limite en dehors desquelles on arrête la simulation (un peu plus grandes que le terrain)

        calcule_acceleration(&a,&forces,*v,cst);
        calcule_position_et_vitesse(x,v,a,*settings);

        //On vérifie tout d'abords si la balle rebondis (z==0?) (on fait ça avant de vérifier si elle est dans le filet car ça peut arriver que la balle rebondisse avant)
        if (x->z<0){ 
            //Réduit les vitesse et inverse la vitesse verticale en supposant que l'angle d'incidence est le même que celui du rebond.
            //effet estéthique plutôt que physique car selon la rotation de la balle elle ne rebondirais pas de cette manière. On considère que les balles ont été frappées plate pour avoir ce genre de rebonds alors que ça n'est pas forcément le cas si on a de la rotation.
            x->z=0;
            v->x*=cst.e; //coeff de restitution des vitesse
            v->y*=cst.e;
            v->z*=-cst.e;
            rebonds+=1;

            //On vérifie si la balle rebondis après avoir dépassé le filet
            if (rebonds == 1 && x->y > 0){ 
                //Annonce dans le terminal des coordonées du point d'impact (si le filet à été dépassé)
                printf("Sans considérer d'aléa, la balle vole pendant %.2f secondes pour atterir en : \n(%.2f,%.2f,0.00)\n",t,x->x,x->y);

                //relevé du point d'impact pour les calculs de probabilités qui suivront si l'option est activée
                impact->x=x->x;
                impact->y=x->y;
                impact->z=x->z;
                //On vérifie si notre coup est IN, le cas échéant on l'annonce out et on annule la partie du calcul faisabilité du service si elle était activée
                if (x->x > limites.x_max || x->x <limites.x_min ||x->y > limites.y_max || x->y < limites.y_min) {
                    printf("Le %s est out.\n",limites.type_de_coup);
                    if (strcmp(settings->calcule_faisabilite,"ON")==0){
                        printf("Nous ne procéderons donc pas au calcul de la faisabilité du %s\n",limites.type_de_coup);
                        settings->calcule_faisabilite="OFF";
                    }
                }
            } 
            //On vérifie si l'on est pas dans le cas particulier ou la balle rebondis avant meme de dépasser le filet (la balle se déplace tjr dans la direction du filet donc v->y >0), si elle revient vers le joueur c'est que y'a eu filet donc on ne peut pas dire rebonds avant filet!
            if (rebonds == 1 && x->y < 0 && v->y>0){
                printf("Votre %s rebondit avant même avoir dépassé le filet!\n",limites.type_de_coup);
                if (strcmp(settings->calcule_faisabilite,"ON")==0) {
                    printf("Le calcul de la faisabilité de ce %s est donc désactivé.\n",limites.type_de_coup);
                    settings->calcule_faisabilite="OFF";
                }
            }
            //rajoute les coordonées du point d'impact au fichier et recommence en haut de la boucle
            fprintf(fichier,"%.5f,%.5f,0.00\n",x->x,x->y);
            t+=settings->dt;
            continue;
        } 

        //Dans le filet? 
        if (filet(intervale_filet,x) == 1) { 
            //On annonce que le coup rentre dans le filet seulement s'il n'a pas déjà rebondis avant de foncer dans celui ci.
            if (rebonds == 0){
                printf("Filet !\n");
                if (strcmp(settings->calcule_faisabilite,"ON")==0){
                    settings->calcule_faisabilite="OFF";
                    printf("Nous ne procéderons donc pas au calcul de la faisabilité du %s.\n",limites.type_de_coup);
                } 
            }
            //On inverse les vitesse et les réduits grandement pour un effet estétique dans le plot de la trajetoire (python)
            v->x=0.1*v->x;
            v->y=-0.1*v->y;
            v->z=0.1*v->z;
            //On suppose que les balles ne tournent plus après le choc
            cst.wx=0;
            cst.wy=0;
            cst.wz=0;
            //On fais ça pour éviter d'entrer dans une boucle infinie car la vitesse devient très faible et on peut se retrouver bloqués dans l'interval filet et donc dans ce if sinon :  
            x->y=-intervale_filet/2;       
            //on ajoute au fichier nos coordonées avant de reprendre en haut de la boucle
            fprintf(fichier,"%.5f,0.00,%.5f\n",x->x,x->z); 
            continue;
        }

        t+=settings->dt;
        //Ajoute les valeurs non impactée par le filet ou un rebonds dans le fichier et recommence la boucle
        fprintf(fichier,"%.5f,%.5f,%.5f\n",x->x,x->y,x->z);
    }

    //fin des calculs des coords de la trajectoire, on ferme le fichier
    fclose(fichier);


    //On vérifie que l'on ne soit pas dans le cas particulier ou l'on sort des limites du terrains sans avoir rebondi. Si c'est le cas, on l'annonce.
    if (rebonds==0){
        printf("Le %s réalisé sort des limites de la simulation sans même attérir!\n",limites.type_de_coup);
        if (strcmp(settings->calcule_faisabilite,"ON")==0){
            settings->calcule_faisabilite="OFF";
            printf("Nous ne procéderons donc pas au calcul de la faisabilité de ce %s étant donné qu'il est largement out.\n",limites.type_de_coup);
        }
    }
}

double gaussian_distributed_number(double randomdomain){
    //Box-Muller transformation (https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
    //rand genere deux valeurs indépendantes suivant la loi uniforme U(0,1) 
    double u1=rand()/randomdomain;
    double u2=rand()/randomdomain;
    //box muller transform : le rayon et l'angle teta
    double radius=sqrt(-2*log(u1));
    double teta=2*M_PI*u2;
    //Calcule and retourne z qui suit la loi normale centrée réduite N(0,1)
    return radius*cos(teta);
}

void calcule_acceleration_avec_alea(struct vecteur * a, struct vecteur * forces, struct vecteur v, struct vecteur vent, struct CST cst){
    //Cette fonction similaire a la fonction calcule_acceleration prends en compte le vent en modifiant la force de frottement de l'air.
    //Plutôt qu'utiliser la vitesse de la balle pour calculer les frottements, on calcule la vitesse relative de celle ci par rapport au vent. On a donc supposé que le vent impact seulement les frottements de l'air.
    double module_vitesse_relative = sqrt((v.x-vent.x)*(v.x-vent.x)+(v.y-vent.y)*(v.y-vent.y)+v.z*v.z); //vent.z = 0
    //si le vent est plus rapide que la balle, alors (v.x-vent.x) < 0 et cf est tjr neg donc la force est positive,ce qui veut dire que la résultante va aller dans la direction de la balle et l'accélerer. Si le v.x-vent.x > 0, alors la force de frottement sera opposée au mvt et moins ou plus importante selon direction du vent.
    //Calcul des forces selons les 3 axes
    forces->x = cst.cf*module_vitesse_relative*(v.x-vent.x)  +  cst.s*(cst.wy*v.z-cst.wz*v.y);     // force de trainee et force de magnus (produit vectoriel)
    forces->y = cst.cf*module_vitesse_relative*(v.y-vent.y)  +  cst.s*(cst.wz*v.x-cst.wx*v.z);     // force de trainee et force de magnus (produit vectoriel)
    
    //La vitesse du vent en z vaut 0 donc pas besoin de modif la force de frott en z car relativement au vent, v_relative_z=v.z
    forces->z = cst.cf*module_vitesse_relative*v.z  +  cst.s*(cst.wx*v.y-cst.wy*v.x) + cst.poids;  // Poids, force trainee et force magnus

    //Deuxième loi newton : somme forces/masse = acceleration
    a->x=forces->x/cst.masse;
    a->y=forces->y/cst.masse;
    a->z=forces->z/cst.masse;
}

struct vecteur simulation_avec_alea(struct vecteur x, struct vecteur  v, struct vecteur a, struct vecteur forces,struct CST cst,struct settings settings, struct vecteur vent) {
    //Version plus simple qu'avant car ne met rien dans fichier et ne s'interesse pas à ce qu'il se passe après rebonds/filet.
    //Elle prends en argument les conditions initiales sur lesquelles un aléa s'est posé et prends en compte le vent également à travers la fonction calcul_acceleration_avec_alea
    struct vecteur impact;
    double t = 0;
    double intervale_filet = v.y*settings.dt;   
    //On lance simulation qui continue du moment ou on reste dans un perimètre acceptable autours du terrain.
    while (x.x <7 && x.x>-7 && x.y< 18 && x.y>-14) { 
        calcule_acceleration_avec_alea(&a,&forces,v,vent,cst);
        calcule_position_et_vitesse(&x,&v,a,settings);
        //On vérifie si l'on a fait filet ; si c'est le cas on retourne l'impact avec -1.
        if (filet(intervale_filet,&x) == 1) {
            impact.x=x.x;
            impact.y=x.y;
            impact.z=-1; //on met -1 pour z pour pouvoir trier les pts après dans fonction de catégorisation des impacts
            return impact; 
        }
        //On vérifie si l'on fait rebonds et si oui on retourne le point d'impact
        if (x.z<0){ 
            impact.x=x.x;
            impact.y=x.y;
            impact.z=0;
            return impact;
        t+=settings.dt;
        }
    }
    //Dans cas rare ou l'on sort du périmètre limite autour du terrains sans rebondire à cause des aléas, on associe la valeur de l'impact avec des -1 partout
    impact.x=-1;
    impact.y=-1;
    impact.z=-1;
    return impact;
}

void calculs_points_impacts_avec_alea(struct vecteur * impacts,struct settings settings,double randomdomain,double module_vitesse_initiale, double teta, double phi,struct alea alea,struct CST * cst, struct vecteur xinitial, struct vecteur vent){
    //Fonction organisatrice de la partie probabilité, elle génére les conditions initiales avec aléas et lance la simulation n fois, tout en stockant chauq epoints d'impact dans le fichier impacts.csv
    
    //On ouvre le fichier dans lequels on dépose tout les points d'impacts
    FILE * fichier = fopen(settings.nom_fichier_impacts, "w");
    fprintf(fichier,"X,Y,Z\n");
    //Le premier point d'impact est celui de la simulation sans aléas, celui ci sera marqué en rouge dans le plot python
    fprintf(fichier,"%.5f,%.5f,%.5f\n",impacts[0].x,impacts[0].y,settings.rayon_de_similitude);

    //il faut stocker les CI de wx,wy et wz car elle seront constamment recalculée et remplacée avec aléas dans la boucle for.
    double wx = cst->wx;
    double wy = cst->wy;
    double wz = cst->wz;

    //on répète la procédure de calcul avec conditions initiales modifiées n fois
    for (int i = 1; i<settings.faisabilite_n+1;i++){
        //On calcule des conditions initiales (variables aléatoires) qui suivent une loi normale centrée en la condition initiale et dont l'écart type est celui rentré dans les options
        //On utilise que N(mu,sigma^2)=sigma*N(0,1)+mu pour génerer nombre aléatoire avec comme espérance la CI et l'écart type ce qu'on a rentré
        //les aléas impactent : les vitesses initiales (angulaires et linéaires) ainsi que les angles donnant la direction du vecteur.
        double v_alea= alea.ecart_type_vitesse*gaussian_distributed_number(randomdomain)+module_vitesse_initiale;
        double teta_alea=alea.ecart_type_angles*gaussian_distributed_number(randomdomain)+teta;
        double phi_alea=alea.ecart_type_angles*gaussian_distributed_number(randomdomain)+phi;
        cst->wx=alea.ecart_type_w*gaussian_distributed_number(randomdomain)+wx;
        cst->wy=alea.ecart_type_w*gaussian_distributed_number(randomdomain)+wy;
        cst->wz=alea.ecart_type_w*gaussian_distributed_number(randomdomain)+wz;

        //On écrase les anciennes positions (de la simulation k-1), vitesse etc pour relancer une nouvelle simulation k , k = 1,...,n
        struct vecteur x = xinitial;                  
        struct vecteur v = coordvitessecartesienne(v_alea,teta_alea,phi_alea);
        struct vecteur a;                         
        struct vecteur forces;

        //On calculera n points d'impacts et les sauveras dans un fichier
        impacts[i] = simulation_avec_alea(x,v,a,forces,*cst,settings,vent);
        fprintf(fichier,"%.5f,%.5f,%.5f\n",impacts[i].x,impacts[i].y,impacts[i].z);
    }
    fclose(fichier);
}

void categorisation_point_impact_avec_alea(struct vecteur* impacts,struct limites limites,struct alea* alea,struct settings settings){
    //Fonction qui catégorise les n impacts des coups avec aléas simulés au préalable et déposé dans le fichier "impacts". 
    alea->dans_cercle=0;
    alea->dans_limites=0;
    alea->dans_limites_et_dans_cercle=0;

    //On catégorise nos n coups affectés par imprécision (on commence a 1 car le 0 est le coup sans aléas)
    for (int i = 1; i<settings.faisabilite_n+1;i++){
        
        if (impacts[i].x >limites.x_min && impacts[i].x < limites.x_max && impacts[i].y > limites.y_min && impacts[i].y < limites.y_max && impacts[i].z != -1){
            //on dénombre les coups dans le terrain/carré de service selon type de coup effectué
            alea->dans_limites+=1;

            //On dénombre les coups dans le terrain et dans le cercle de similitude
            if (sqrt(pow((impacts[i].x-impacts[0].x),2) + pow((impacts[i].y-impacts[0].y),2)) < settings.rayon_de_similitude && impacts[i].z != -1){
                alea->dans_limites_et_dans_cercle+=1;
            }
        }
        //on dénombre les coups dans le cercle de similitude (et pas dans le filet)
        if (sqrt(pow((impacts[i].x-impacts[0].x),2) + pow((impacts[i].y-impacts[0].y),2)) < settings.rayon_de_similitude && impacts[i].z != -1){
            alea->dans_cercle+=1;
        }
    }
}

//MAIN

int main(int argc, char * argv[]) {
    struct CST cst;
    struct vecteur xinitial; //vecteur condition initales
    struct settings settings;
    struct alea alea;

    //paramètres par défaut (même que sur site internet), au cas ou on ne complète pas tout sur site

    xinitial.x =0.5                                          ;//[m]                      voir rendu pour placement axes
    xinitial.y =-11.89                                       ;//[m]                      Ligne du fond à -11.89 m 
    xinitial.z = 2.5                                         ;//[m]
    double module_vitesse_initiale = 120                 /3.6;//[km/h]->[m/s]
    double teta = 90                                *M_PI/180;//[degré]->[radian]        entre 0 et 90 vers haut, entre 90 et 180 vers bas.
    double phi  = 100                               *M_PI/180;//[degré]->[radian]        entre 0 et 90 vers la droite,90 et 180 vers la gauche.
    //vitesse angulaire      
    cst.wx = -1500                                 *2*M_PI/60;//[rpm (revolution per minute)]->[rad/s]     
    cst.wz = 400                                   *2*M_PI/60;      
    cst.wy = -400                                  *2*M_PI/60;      
  
    //Aléas  
    alea.ecart_type_vitesse=5                            /3.6;//[km/h]->[m/s]        La précision du joueurs sur les conditions initiales 
    alea.ecart_type_angles=1.5                      *M_PI/180;//[degré]->[radian] 
    alea.ecart_type_w=75                           *2*M_PI/60;//[rpm]->[rad/s]
    alea.vent= 6                                         /3.6;//[km/h]->[m/s]        
    alea.phi_vent=45                                *M_PI/180;//[degré]->[radian]    angle entre vecteur vitesse du vent et axe x dans notre système de coordonée

    //Settings
    settings.dt=0.001                                       ;//incrément de la simulation (un dt trop petit rends la simulation plus précise mais plus longue) ; réglage initial conseillé : 0.001
    settings.type_de_coup ="service_gauche"                 ;//"service_droite","service_gauche" ou "coup_normal"   || important pour déterminer les limites dans lesquels la balle doit tomber pour ne pas être out. si on donne pas une des trois options possible, on est sur un "coup normal" par défaut.
    settings.frottement_air = "ON"                          ;//Prendre en compte les frottements de l'air 'ON' ou 'OFF', si autre chose on a du ON par défaut
    //ATTENTION : On considère que le vent n'impacte que la force de frottement de l'air, donc si on la desactive dans les settings le vent n'aura pas d'effet non plus.
    
    settings.calcule_faisabilite ="ON"                 ;//Si l'on veut calculer la faisabilité du service rentré avec les conditions initiales (basé sur notre modèle stochastique) ; 'ON' ou 'OFF' et si autre chose 'OFF par défaut
    settings.faisabilite_n=50000                       ;//Nombre de fois qu'on répètera simulations avec aléas pour faire des calculs de probabilité ensuite (50000 recommandé), le plus grand n est, plus l'estimation de la probabilité converge vers la prob théorique
    settings.rayon_de_similitude=0.75                  ;//[m] Rayon du cercle centré en le point d'impact du coup sans aléa dans lequel ont considère que les coups qu'on a efféctués avec aléa sont similaires au coups sans aléa. (on considère pas les filets)

    settings.nom_fichier = "trajectoire.csv"           ;//Nom du fichier dans lequel les coordonées de la trajectoire sans aléas sont entrées. extension csv important
    settings.nom_fichier_impacts="impacts.csv"         ;//Nom du fichier dans lequel on pose nos n point d'impact simulés avec aléas. extension csv important


    // Lire les paramètres du site web

	char * parametres = argc > 1 ? argv[1] : NULL;
	while (parametres != NULL) {
		char * equal = strchr(parametres, '=');
		if (equal == NULL) break;
		equal[0] = 0;
		char * ampersand = strchr(equal + 1, '&');
		if (ampersand != NULL) ampersand[0] = 0;
		char * key = parametres;
		char * value = equal + 1;

        //On converti les paramètre en double/int/string et les convertis en les bonnes unités (par ex km/h -> m/s)

		if (strcmp(key, "x") == 0) {
			xinitial.x = atof(value);
		} else if (strcmp(key, "y") == 0) {
			xinitial.y = atof(value);
		} else if (strcmp(key, "z") == 0) {
			xinitial.z = atof(value);
		} else if (strcmp(key, "v") == 0) {
			module_vitesse_initiale = atof(value)/3.6;
		} else if (strcmp(key, "teta") == 0) {
			teta = atof(value)*M_PI/180;
		} else if (strcmp(key, "phi") == 0) {
			phi = atof(value)*M_PI/180;
		} else if (strcmp(key, "wx") == 0) {
			cst.wx = atof(value)*2*M_PI/60;
		} else if (strcmp(key, "wy") == 0) {
			cst.wy = atof(value)*2*M_PI/60;
		} else if (strcmp(key, "wz") == 0) {
			cst.wz = atof(value)*2*M_PI/60;
		} else if (strcmp(key, "ec_v") == 0) {
		    alea.ecart_type_vitesse = atof(value)/3.6;
		} else if (strcmp(key, "ec_a") == 0) {
			alea.ecart_type_angles = atof(value)*M_PI/180;
		} else if (strcmp(key, "ec_w") == 0) {
			alea.ecart_type_w = atof(value)*2*M_PI/60;
		} else if (strcmp(key, "vent") == 0) {
			alea.vent = atof(value)/3.6;
		} else if (strcmp(key, "phi_vent") == 0) {
			alea.phi_vent = atof(value)*M_PI/180;
		} else if (strcmp(key, "type") == 0) {
			settings.type_de_coup = value;
		} else if (strcmp(key, "frott") == 0) {
			settings.frottement_air = value;
		} else if (strcmp(key, "proba") == 0) {
			settings.calcule_faisabilite = value;
		} else if (strcmp(key, "n") == 0) {
			settings.faisabilite_n = atoi(value);
		} else if (strcmp(key, "r") == 0) {
			settings.rayon_de_similitude = atof(value);
		}
		parametres = ampersand == NULL ? NULL : ampersand + 1;
	}

    //Définition des constantes spécifiques au tennis et dernier settups automatiques indispensables pour utiliser les fonctions
    cst.masse = 0.057;                             //[KG] 
    cst.poids = cst.masse*-9.81;                   //[N], selon axe z
    cst.rho = 1.2;                                 //[kg/m^3]
    cst.r = 0.0326;                                //[m]
    cst.s = 0.5*M_PI*cst.rho*cst.r*cst.r*cst.r;    //[kg]   il s'agit du facteur de force magnus
    cst.cf = -0.5*cst.rho*M_PI*cst.r*cst.r*0.5;    //[kg/m] coefficient frottement Cx = 1/2 pour air
    cst.e = 0.8;                                   //coefficient de restitution des vitesses après collision avec le sol ( /!\ approx que le coup est plat)
    
    //implémentaton structures nécessaires pour utilisation fonctions
    struct vecteur x = xinitial;                //vecteur position, on conserve le vecteur position initiale pr la partie proba après si elle est activée
    struct vecteur v = coordvitessecartesienne(module_vitesse_initiale,teta,phi); //vecteur vitesse
    struct vecteur a;                           //vecteur accéleration
    struct vecteur forces;                      //vecteur force résultante
    struct vecteur impact_sans_alea;            //coord du point d'impact sans aléa
    struct limites limites = determination_limites(settings.type_de_coup);
    if (strcmp(settings.frottement_air,"OFF")==0) cst.cf = 0;
    //On lance les fonctions depuis MAIN
    
    //Simulation sans alea pour coup principal
    printf("\n--------\nSimulation du %s sans considérer les aléas:\n\n",limites.type_de_coup);
    simulation_sans_alea(&x,&v,a,forces,cst,&settings,&impact_sans_alea,limites);

    //Si le coup est réussi et la faisabilité est "ON", on calcule la faisabilité de celui-ci selon la précision rentrée en conditions initiales
    
    
    if (strcmp(settings.calcule_faisabilite,"ON")==0){
        //initialise les fonctions d'ou sortiront nos nombres aléatoires
        int myTime = time(NULL);
        srand(myTime);
        double randomdomain = RAND_MAX+1.0;
        //Le fait de lancer dans le vide deux fois rand fait que l'aléatoire fonctionne mieux ensuite (sous windows en tout cas)
        rand();
        rand();

        //Créer un tableau de structure dont la taille permet de mettre tous les points d'impacts avec aléas qu'on veut
        struct vecteur impacts[settings.faisabilite_n+1];

        //on converti l'aléa du vent en vecteur vitesse (en considérant vent.z = 0,vent dans plan horizontal seulement -> on met teta = pi/2 ds notre fonction de conversion)
        struct vecteur vent = coordvitessecartesienne(alea.vent,M_PI/2,alea.phi_vent);

        //le premier point d'impact est celui sans aléas
        impacts[0]=impact_sans_alea;

        //On remplis ce tableau de structures avec nos n répétitions de simulation avec aléas
        calculs_points_impacts_avec_alea(impacts,settings,randomdomain,module_vitesse_initiale,teta,phi,alea,&cst,xinitial,vent);

        //On catégorise les n coups avec aléas
        categorisation_point_impact_avec_alea(impacts,limites,&alea,settings);

        //On calcule les probas, leur descriptions est donnée dans le terminal
        double P1=(double) alea.dans_cercle/settings.faisabilite_n;
        double P2=(double) alea.dans_limites/settings.faisabilite_n;
        double P3=(double) alea.dans_limites_et_dans_cercle/settings.faisabilite_n;
        double P4=1-P2;

        //On imprime les résultats des calculs de probabilité
        printf("\n\nCalculs empiriques de la faisabilité de ce %s selon la précision \nrentrée en condition initiale (à partir de %d observations):\n\n******\n",limites.type_de_coup,settings.faisabilite_n);
        printf("\nProbabilité de faire un %s qui ratterisse dans le %s (in): %.3f\n",limites.type_de_coup,limites.nom,P2);
        printf("\nProbabilité de faire un %s similaire au %s sans aléa \n(qui atterrisse dans un rayon de %.2f m par rapport au coup sans aléa) : %.3f\n",limites.type_de_coup,limites.type_de_coup,settings.rayon_de_similitude,P1);
        printf("\nProbabilité de faire un %s similaire au %s sans aléas (dans le cercle de similitude) \net dans le %s (in): %.3f\n",limites.type_de_coup,limites.type_de_coup,limites.nom,P3);
        printf("\nProbabilité de faire un filet ou un %s out: %.3f\n",limites.type_de_coup,P4);
        printf("\n******\n\nJetez un oeil sur les graphs ci-dessous illustrant votre coup et votre précision!\n\n");
    }
    else {
        //On imprime ligne de fin de programme si l'on ne passe pas par le calcul de la faisabilité du coup.
        printf("\nLe calcul de la faisibilité étant désactivé,la faisabilité du coup ne sera pas évaluée.\nVous pouvez toutefois retrouver le plot de la trajectoire de votre coup ci-dessous.\n--------");
        //Fait en sorte qui le prog python ne calcule pas le troisième graph.
        FILE * fichier = fopen(settings.nom_fichier_impacts, "w");
        fprintf(fichier,"X,Y,Z\n");
        fprintf(fichier,"0.0,0.0,0.0");
        fclose(fichier);
    }
    //fait tri sur serveur pour éviter de le surcharger d'images non utilisées.
    system("rm image/*");
    system("./graph_versionsite.py");
    return 0;
}