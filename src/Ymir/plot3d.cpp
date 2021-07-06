#include "plot3d.h"
#include <map> //for overlay of heatmaps
/** See plot3d.h for infos:
 *  Displays structures, superProteins or positions in space as spheres.
 *  - use #define ALLOW_GRAPHICS in plot3d.h to use. The library freeglut is then needed. */

#ifdef ALLOW_GRAPHICS
#include "proteins.h"
#include "lattice.h"
#include "compact.h"
#include "zaptrackball.h"
#include "../Tools/zaprandom.h"
#include <sstream>
#include <iostream>
#include <iomanip> // setprecision
#include <fstream>
using namespace std;

void display();
void init();

// IL/SOIL library is used to export glut display as pictures. Beware of capital letters in SOIL.h
#include "soil/SOIL.h"
// IL was an older implementation of SOIL. also possible.
/* #include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>*/

// Storage of all elemets on display (to be redisplay anytime) - given as pointers, should not be destroyed by the user!
static vector<superProtein*> currentProts;
static vector<bool> keepAlways;                     // for each superProt, says if always present or not when using d/q
static vector<set<int>*> currentSurfaces;
static vector<vector<double>> currentPoints;
static vector< vector< std::pair<int, double> >* > currentHeatmap; // map: position => something between 0 and 1
static size_t IDheatmap = 0;

// Variables to represent the camera/view
//int spinning = 0,
static bool moving = 0;
static int beginx, beginy;     // saves beginning position when clicked move
static int W = 300, H = 300; // resizing needs to be done ...
static double moveX = 0, moveY = 0, moveZ = 0;
static double manualPointX = XWidth / 2, manualPointY = YWidth / 2, manualPointZ = ZWidth / 2;

static float curquat[4];   // current quaternion
static float lastquat[4];

static GLdouble bodyWidth = 3.0;
static bool newModel = 1;
//static int scaling = 4; // 1
static double scalefactor = 5; //1.0;
static bool waitRightClick = false;
static int indexStructDisplayed = -1;
static double sizeSphere = 0.75; //0.2;
static bool showAAColors = true; // false
static bool showHeatmap = false;
static bool visualizeInterface = false;
//static vector<double> heatmapZeroColor = {0.1, 0.1, 0.1};
//static vector<double> heatmapOneColor = {1.0, 0.7, 0.7};
//static vector<double> heatmapCoreColor = {35./256., 220./256., 170./256.};
static vector<double> heatmapZeroColor = {0.1, 0.1, 0.1};
static vector<double> heatmapOneColor = {149./256., 196./256., 234./256.};
static vector<double> heatmapCoreColor = {149./256., 196./256., 234./256.};
//static vector<double> defaultProteinColor = {220./256., 38./256., 127./256.}; // Magenta
// static vector<double> defaultProteinColor = {254./256., 97./256., 0./256.}; // Orange
//static vector<double> defaultProteinColor = {255./256., 176./256., 0./256.}; // Gold
//static vector<double> defaultProteinColor = {241./256., 41./256., 138./256.}; // 1FBI
//static vector<double> defaultProteinColor = {201./256., 148./256., 199./256.};  //1ADQ
//static vector<double> defaultProteinColor = {206./256., 18./256., 86./256.};  //1H0D
static vector<double> defaultProteinColor = {256./256., 256./256., 256./256.};  //1H0D


static bool showAxis = false;
static bool showForbPositions = true;
static bool saving = false;
static bool backgroundColorBlack = false;
static string PictureFileNamePrefix = "";

void setPictureFileNamePrefix(string newPrefix){
    PictureFileNamePrefix = newPrefix;
}

void print3DState(){
    cout << /*"Spinning " << spinning <<*/ ", moving " << moving << ", beginx " << beginx << ", beginy " << beginy << ", bodyWidth " << bodyWidth << ", scaleFactor " << scalefactor << endl;
}

void waitForRightClick(){
    waitRightClick = true;
    while(waitRightClick){
        glutMainLoop();
    }
}

void recalcModelView(void){
    GLfloat m[4][4];
    glPopMatrix();    //???
    glPushMatrix();
    build_rotmatrix(m, curquat);
    glMultMatrixf(&m[0][0]);
    if (scalefactor == 1.0) {
        glDisable(GL_NORMALIZE);
    } else {
        glEnable(GL_NORMALIZE);
    }
    glScaled(scalefactor, scalefactor, scalefactor);
    //glTranslatef(-8, -8, -bodyWidth / 2);
    newModel = 0;
}

void myReshape(int w, int h){
    glViewport(0, 0, w, h);
    W = w;
    H = h;
}

void mouse(int button, int state, int x, int y){
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        //spinning = 0;
        glutIdleFunc(nullptr);
        moving = 1;
        beginx = x;
        beginy = y;
        // glutPostRedisplay will be catched later because of motion function !!
        /*if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
          scaling = 1;
        } else {
          scaling = 0;
        }*/
    }
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        moving = 0;
    }
    if(button == 3){
        //cout << "scroll up" << endl; // note: always called twice
        scalefactor *= 1.08;
        newModel = 1;             // needs to recompute the model (object)
        glutPostRedisplay(); // needs
    }
    if(button == 4){
        //cout << "scroll down" << endl;
        scalefactor *= 0.92;
        newModel = 1;
        glutPostRedisplay();
    }
}

// void animate(void){
//  add_quats(lastquat, curquat, curquat);
//  newModel = 1;
//  glutPostRedisplay();
//}

void motion(int x, int y)
{
  /*if (scaling) {
    //scalefactor = scalefactor * (1.0 + (((float) (beginy - y)) / H));
    beginx = x;
    beginy = y;
    newModel = 1;
    glutPostRedisplay();
    return;
  }*/
  if (moving) {
    trackball(lastquat,
      (static_cast<float>(2.0) * beginx - W) / W,
      (H - static_cast<float>(2.0) * beginy) / H,
      (static_cast<float>(2.0) * x - W) / W,
      (H - static_cast<float>(2.0) * y) / H
      );
    beginx = x;
    beginy = y;
    //spinning = 1;

    /////glutIdleFunc(animate);
    /// instead, does it only once
    add_quats(lastquat, curquat, curquat);
    newModel = 1;
    glutPostRedisplay();
  }
}

static GLboolean lightZeroSwitch = GL_TRUE, lightOneSwitch = GL_TRUE;

void controlLights(int value){
    switch (value) {
    case 1:
        lightZeroSwitch = !lightZeroSwitch;
        if (lightZeroSwitch) {
            glEnable(GL_LIGHT0);
        } else {
            glDisable(GL_LIGHT0);
        }
        break;
    case 2:
        lightOneSwitch = !lightOneSwitch;
        if (lightOneSwitch) {
            glEnable(GL_LIGHT1);
        } else {
            glDisable(GL_LIGHT1);
        }
        break;
#ifdef GL_MULTISAMPLE_SGIS
    case 3:
        if (glIsEnabled(GL_MULTISAMPLE_SGIS)) {
            glDisable(GL_MULTISAMPLE_SGIS);
        } else {
            glEnable(GL_MULTISAMPLE_SGIS);
        }
        break;
#endif
    case 4:
        glutFullScreen();
        break;
    case 5:
        exit(0);
    }
    glutPostRedisplay();
}

void vis(int visible)
{
  if (visible == GLUT_VISIBLE) {
    //if (spinning)
    //  glutIdleFunc(animate);
  } else {
    //if (spinning)
    //  glutIdleFunc(NULL);
  }
}

void drawBitmapText(string s,float x,float y,float z){
    glRasterPos3f(x, y,z);
    for(unsigned int i = 0; i < s.size(); ++i)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, static_cast<int>(s[i]));
}

void drawBitmapText2D(string s,float x,float y){
    glRasterPos2f(x, y);
    for(unsigned int i = 0; i < s.size(); ++i)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, static_cast<int>(s[i]));
}


void SpecialInput(int key, int x, int y){
    cout << "Special input " << key << " x=" << x << " y=" << y << endl;
    switch(key){
    case GLUT_KEY_UP:
        moveY += 1.;
        break;
    case GLUT_KEY_DOWN:
        moveY -= 1.;
        break;
    case GLUT_KEY_LEFT:
        moveX -= 1.;
        break;
    case GLUT_KEY_RIGHT:
        moveX += 1.;
        break;
    }
display();
}

void keyboard(unsigned char key, int /*x*/, int /*y*/){

    switch(key) {
        case 'q': case 'Q':
            indexStructDisplayed--;
            if(indexStructDisplayed < -1) indexStructDisplayed = static_cast<int>(keepAlways.size()) - 1;
        break;
        case 'd': case 'D':
            indexStructDisplayed++;
            //cout << keepAlways.size() << endl;
            if(indexStructDisplayed >= static_cast<int>( keepAlways.size())) indexStructDisplayed = -1;
        break;
        case 'l':
            sizeSphere /= 1.2;
        break;
        case 'L':
            sizeSphere *= 1.2;
        break;
        case 'A': case 'a':
            showAAColors = !showAAColors;
            if(showAAColors) showHeatmap = false;
        break;
        case 'H': case 'h':
            showHeatmap = !showHeatmap;
            if(showHeatmap) showAAColors = false;

        break;
        case 'J': case 'j':
            IDheatmap++;
            if(IDheatmap >= currentHeatmap.size() + 1) IDheatmap = 0;
        break;
        case 'B': case 'b':
            backgroundColorBlack = !backgroundColorBlack;
            if(backgroundColorBlack){
                glClearColor(0.0, 0.0, 0.0, 1.0);
            } else {
                glClearColor(1.0, 1.0, 1.0, 1.0);
            }
        break;
        case 'S': case 's':
            showAxis = !showAxis;
        break;
        case 'F': case 'f':
            showForbPositions = !showForbPositions;
        break;

        case 'X':
            manualPointX += 1;
        break;
        case 'Y':
            manualPointY += 1;
        break;
        case 'Z':
            manualPointZ += 1;
        break;
        case 'x':
            manualPointX -= 1;
        break;
        case 'y':
            manualPointY -= 1;
        break;
        case 'z':
            manualPointZ -= 1;
        break;
        case 'P': case 'p':
            cout << "(" << manualPointX << "," << manualPointY << "," << manualPointZ << ")" << lattice::idFromPosisition(manualPointX, manualPointY, manualPointZ) << endl;
        break;
        case 'O': case 'o':
            saving = !saving;
            if(saving) cout << "... Saving images every time refreshed " << endl;
            else cout << "... Stopping saving images" << endl;
        break;

        case 'I': case 'i':
            visualizeInterface = !visualizeInterface;
        break;
    case '+':
        scalefactor *= 1.08;
        newModel = 1;             // needs to recompute the model (object)
        glutPostRedisplay(); // needs
        break;
    case '-':
        scalefactor /= 1.08;
        newModel = 1;             // needs to recompute the model (object)
        glutPostRedisplay(); // needs
        break;
        default:
        break;
    }
    glutPostRedisplay(); /* this redraws the scene without
        waiting for the display callback so that any changes appear
        instantly */
}

// This is actually the main init function, to be called.
void glDisplay(int argc, char** argv, bool append){
    // The content of this function should only be called once
    static bool alreadyInit = false;
    if(alreadyInit) {
        if(!append) clearDisplay();
        return;
    }
    alreadyInit = true;

    // creating argc and argv if not given.
    if(!argv){
        //char* a[] = {(char*) "banana", (char*) "apple"}; // needs an arbitrary argc and argv
        char* a = (char*) "a";
        char* v[] = {a};
        argv = v;
    }

    //Initialize windows
    glutInit(&argc, argv);

    //Single = one buffer ; double = work on a hidden buffer and then swap
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutInitWindowSize(700, 700);
    glutInitWindowPosition(200,200);

    // Initializes the rotation matrices
    trackball(curquat, 0, 0, 0, 0);
    //trackball(curquat, 0.3, -0.4, 0.7, 0.3);

    // Name of the window
    glutCreateWindow("3D structures showbox");

    // Defines the functions that need to be called when user does things
    glutDisplayFunc(display);       // whenever needs to redisplay
    glutReshapeFunc(myReshape);     // when size of window is changed
//  glutVisibilityFunc(vis);        // ???
    glutMouseFunc(mouse);           // When mouse does something
    // glutMouseWheelFunction(mouse); doesn't exist on the official glut => get button 3 and 4
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(SpecialInput);  // for arrows

    // If want to have menus.
    /* glutCreateMenu(controlLights);
    glutAddMenuEntry("Toggle right light", 1);
    glutAddMenuEntry("Toggle left light", 2);
    glutAddMenuEntry("Full screen", 3);
    glutAddMenuEntry("Quit", 4);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    glEnable(GL_CULL_FACE); */

    // in case you want to use Devil library
    /*ilInit();
    iluInit();
    ilutRenderer(ILUT_OPENGL);*/

    // Function (see below) to setup the object and the camera. Can be called multiple times
    init();

    // the display function is called automatically if sth needs to be redisplayed.
    //glutDisplayFunc(display);//(somatic, CTL));
    //glutMainLoop();
}


void init(void)
{
    //Background color. Last 0.0 is for transparency
    glClearColor(1.0, 1.0, 1.0, 1.0);
    // that would be white; glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // to specify another position for the light
    GLfloat lightOnePosition[] = {0, 0, ZWidth, 0.0};
    glLightfv(GL_LIGHT0, GL_POSITION, lightOnePosition);
    GLfloat lightOneColor[] = {1.0, 1.0, 1.0, 1.0};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightOneColor);    // very important
    // OR/AND
    //GLfloat specular[] = {1.0f, 1.0f, 1.0f , 1.0f};
    //glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    // OR/AND
    //GLfloat ambient[] = {1.0, 1.0, 1.0, 1.0};
    //glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

    #define zoom 1
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, 1.0, 1, 500); // Note : deph test works only if the first plane is > 0
    glMatrixMode(GL_MODELVIEW);
    gluLookAt(0, 0, ZWidth * 2.0 /*+ ZWidth*/,  /* eye is at (0,0,30) */
      0, 0, 0,      /* center is at (0,0,0) */
      0.0, 1.0, 0.);      /* up is in positive Y direction */
    //glOrtho(-0.1, 0.1, -0.1, 0.1, -0.1, 0.1);
    //glPushMatrix(); /// Phi: might need to remove
    glEnable(GL_COLOR_MATERIAL); // very important, if not there is no color at all
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glEnable(GL_DEPTH_TEST);
}






void addToDisplay(superProtein* new_s, bool always_visible){
    if(!new_s) {cerr << "ERR: addToDisplay(NULL)" << endl; return;}
    keepAlways.push_back(always_visible);
    //for(unsigned int i = 0; i < currentProts.size(); ++i){
    //    if(currentProts[i] == new_s) return;
    //}
    currentProts.push_back(new_s);
}

void addToDisplay(struct3D* new_s, bool always_visible){
    superProtein* protForDisplay = new superProtein(*new_s);
    addToDisplay(protForDisplay, always_visible);
}

void addDissected(superProtein* new_s, moveDirection startDir, int startPos){
    // Note: done from central position and undefined start.
    if(!new_s) {cerr << "ERR: addDissected(NULL)" << endl; return;}
    if(!new_s->structure) {cerr << "ERR: addDissected(NULL structure)" << endl; return;}
    string S = new_s->structure->sequence;
    size_t L = S.size();
    // Less than 2 moves doesn't need dissect
    if(L <= 2){
        addToDisplay(new_s, false);
    }
    for(size_t i = 0; i < L; ++i){ //££tag
        superProtein* s2 = new superProtein(new_s->structure->sequence.substr(0,i+1), startPos);
        addToDisplay(s2, false);
    }
}



void addToDisplay(set<int>* new_s){
    for(unsigned int i = 0; i < currentSurfaces.size(); ++i){
        if(currentSurfaces[i] == new_s) return;
    }
    currentSurfaces.push_back(new_s);
}

void addToDisplay(vector<vector<double>>* listpoints){
    if(!listpoints) return;
    cerr << "Clear POINTS" << endl;
    currentPoints.clear();
    for(unsigned int i = 0; i < listpoints->size(); ++i){
        currentPoints.push_back((*listpoints)[i]);
    }
}

// this is not consistent...
void addToDisplay(vector< std::pair<int, double> > * newHeatmap){
    currentHeatmap.push_back(newHeatmap);
}

bool isin(double x, double vmin, double vmax){
    return((x >= vmin) && (x <= vmax));
}

void setHeatmapColorZero(vector<double> v){
    if(v.size() != 3){cerr << "ERR: setHeatmapColorForZero should be vector of size 3" << endl; return;}
    if((isin(v[0], 0, 1)) && (isin(v[1], 0, 1)) && (isin(v[2], 0, 1))){
        heatmapZeroColor = v;
    } else {
        cerr << "ERR: setHeatmapColorZero, The colors for the heatmap should lie between 0 and 1" << endl;
    }
}

void setHeatmapColorOne(vector<double> v){
    if(v.size() != 3){cerr << "ERR: setHeatmapColorOne should be vector of size 3" << endl; return;}
    if((isin(v[0], 0, 1)) && (isin(v[1], 0, 1)) && (isin(v[2], 0, 1))){
        heatmapOneColor = v;
    } else {
        cerr << "ERR: setHeatmapColorOne, The colors for the heatmap should lie between 0 and 1" << endl;
    }
}


void clearDisplay(){
    //cerr << "WHY cleared" << endl;
    currentProts.clear();
    currentSurfaces.clear();
    // Philippe: this needs to be restored !!! - heatmaps need to be cleared as well
    //currentPoints.clear();
    currentHeatmap.clear();
}


void display(){
    //cout << "Display - " << indexStructDisplayed << endl;
    if (newModel)
      recalcModelView();

    // Question : how to write text on the screen ???
    /*glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glColor3f(1.0f, 0.0f, 0.0f);//needs to be called before RasterPos
    text ???
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();*/

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    //when drawing something, give color first, then draw
    //f signifies float
    glLineWidth(3.0);
    glColor4d(0.45, 0.0, 0.45, 1.0);
    glColor4d(1, 1, 1, 1.0);

    glTranslated(moveX, moveY, moveZ);

    if(showAxis){
        // Borders of the cube
        for(int ix = 0; ix <=1; ++ix){
            for(int iy = 0; iy <=1; ++iy){
                for(int iz = 0; iz <=1; ++iz){
                    glBegin(GL_LINES);
                    glVertex3d(XWidth * (ix-0.5) , YWidth * (iy-0.5) , -ZWidth * (iz-0.5) );
                    glVertex3d(XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
                    glVertex3d(XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
                    glVertex3d(XWidth * (ix-0.5) , - YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
                    glVertex3d(XWidth * (ix-0.5) , - YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
                    glVertex3d(XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
                    glVertex3d(XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
                    glVertex3d(- XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
                    glEnd();
                }
            }
        }

        // Axes
        glColor4d(0.3, 0.3, 0.85, 1.0);
        for(int i = 0; i < XWidth ; ++i){
            glBegin(GL_LINES);
            glVertex3d(i -0.1 - 0.5*XWidth, 0 , 0 );
            glVertex3d(i + 0.1 - 0.5*XWidth, 0 , 0 );
            glEnd();
            glBegin(GL_LINES);
            glVertex3d(0,i - 0.1 - 0.5*YWidth, 0  );
            glVertex3d(0,i + 0.1 - 0.5*YWidth, 0 );
            glEnd();
            glBegin(GL_LINES);
            glVertex3d(0,0,i - 0.1  - 0.5*ZWidth);
            glVertex3d(0,0,i + 0.1 - 0.5*ZWidth);
            glEnd();
        }
        drawBitmapText("X",XWidth / 2.0,0,0);
        drawBitmapText("Y",0,YWidth / 2.0,0);
        drawBitmapText("Z",0,0,ZWidth/2.0);
    }

    #define randomCoeff 0
    //#define sizeSphere 0.2
    for(size_t is = 0; is < currentProts.size(); ++is){
        if((indexStructDisplayed == -1) || (indexStructDisplayed == static_cast<int> (is)) || (keepAlways[is])){
            superProtein* s = currentProts[is];
            if(s == nullptr) continue;
            stringstream tt;
            tt << "S="<< s->structure->sequence;
            drawBitmapText(tt.str().c_str(),-XWidth / 3.0, -YWidth/2.0 -is*5 , 0);

            size_t nR = s->points.size();
            vector<int> posPrec;
            vector<int> pos;

            if(keepAlways[is]){                      
//              glColor4d(1.0, 1.0, 1.0, 1.0);
                glColor4d(defaultProteinColor[0], defaultProteinColor[1], defaultProteinColor[2], 1.0);
            } else {
                double dr = random::uniformDouble(0,1);
                double dg = random::uniformDouble(0,1);
                double db = random::uniformDouble(0,1);
                glColor4d(dr, dg, db, 0.2);
            }

            double dx1 = random::uniformDouble(0,randomCoeff*0.2);
            double dy1 = random::uniformDouble(0,randomCoeff*0.2);
            double dz1 = random::uniformDouble(0,randomCoeff*0.2);
            for(size_t i = 0; i < nR; ++i){
                posPrec = pos;
                pos = lattice::positionFromID(s->points[i].IDposition);

                if(showAAColors){
                   // Amino acid colors rasmol
//                    switch(s->points[i].TypeResidue){
//                        case Asp: case Glu:             glColor4d(230./255.,230./255., 10./255.,1.0); break;
//                        case Cys: case Met:             glColor4d(230./255.,230./255.,  0./255.,1.0); break;
//                        case Arg: case Lys:             glColor4d( 20./255., 90./255.,255./255.,1.0); break;
//                        case Thr: case Ser:             glColor4d(250./255.,150./255.,  0./255.,1.0); break;
//                        case Tyr: case Phe:             glColor4d( 50./255., 50./255.,170./255.,1.0); break;
//                        case Asn: case Gln:             glColor4d(  0./255.,220./255.,220./255.,1.0); break;
//                        case Gly:                       glColor4d(235./255.,235./255.,235./255.,1.0); break;
//                        case Ile: case Leu: case Val:   glColor4d( 15./255.,130./255., 15./255.,1.0); break;
//                        case Ala :                      glColor4d(200./255.,200./255.,200./255.,1.0); break;
//                        case Trp:                       glColor4d(180./255., 90./255.,180./255.,1.0); break;
//                        case His:                       glColor4d(130./255.,130./255.,210./255.,1.0); break;
//                        case Pro:                       glColor4d(220./255.,150./255.,130./255.,1.0); break;
//                        case UndefinedYet:              glColor4d(0.1, 0.1, 0.1, 1.0); break;
//                        case NB_AAs:                    glColor4d(0.05, 0.05, 0.05, 1.0); break;
//                    }
                    // shapely colors rasmol
                    switch(s->points[i].TypeResidue){
                        case Ala:             glColor4d(140./255.,  255./255.,   140./255.,1.0); break;
                        case Gly:             glColor4d(255./255.,  255./255.,   255./255.,1.0); break;
                        case Leu:             glColor4d( 69./255.,   64./255.,    69./255.,1.0); break;
                        case Ser:             glColor4d(255./255.,  112./255.,    66./255.,1.0); break;
                        case Val:             glColor4d(255./255.,  140./255.,   255./255.,1.0); break;
                        case Thr:             glColor4d(184./255.,   76./255.,     0./255.,1.0); break;
                        case Lys:             glColor4d( 71./255.,   71./255.,   184./255.,1.0); break;
                        case Asp:             glColor4d(160./255.,    0./255.,   166./255.,1.0); break;
                        case Ile:             glColor4d(  0./255.,   76./255.,     0./255.,1.0); break;
                        case Asn:             glColor4d(255./255.,  124./255.,   112./255.,1.0); break;
                        case Glu:             glColor4d(102./255.,    0./255.,     0./255.,1.0); break;
                        case Pro:             glColor4d( 82./255.,   82./255.,    82./255.,1.0); break;
                        case Arg:             glColor4d(  0./255.,    0./255.,   124./255.,1.0); break;
                        case Phe:             glColor4d( 83./255.,   76./255.,    66./255.,1.0); break;
                        case Gln:             glColor4d(255./255.,   76./255.,    76./255.,1.0); break;
                        case Tyr:             glColor4d(140./255.,  112./255.,    76./255.,1.0); break;
                        case His:             glColor4d(112./255.,  112./255.,   255./255.,1.0); break;
                        case Cys:             glColor4d(255./255.,  255./255.,   112./255.,1.0); break;
                        case Met:             glColor4d(184./255.,   160./255.,   66./255.,1.0); break;
                        case Trp:             glColor4d( 79./255.,    70./255.,    0./255.,1.0); break;
                        case UndefinedYet:    glColor4d(0.1, 0.1, 0.1, 1.0); break;
                        case NB_AAs:          glColor4d(0.05, 0.05, 0.05, 1.0); break;
                    }
                }
                //GLfloat d1[] = { 0.1, 0.4, 0.9, 1.0 };
                //glMaterialfv(GL_FRONT,GL_DIFFUSE,d1);
                //glEnable(GL_COLOR_MATERIAL);
                //glColor3f(0.1, 0.4, 0.9);




//                    if( i < 2) // first 2 positions are different
//                        glColor4d(0.8, 1.0, 0.8, 1.0);
//                    else
//                        glColor4d(1.0, 0.3, 0.3, 1.0);

                glTranslated(dx1 + pos[0] - XWidth / 2., dy1 + pos[1] -YWidth / 2., dz1 + pos[2] - ZWidth / 2.);
                glutSolidSphere(sizeSphere, 50, 50);
                glTranslated(-dx1 -pos[0] + XWidth / 2., -dy1 -pos[1] +YWidth / 2., -dz1 -pos[2] + ZWidth / 2.);
                if(i > 0){
                    glBegin(GL_LINES);
                    glVertex3d(dx1 + posPrec[0] - XWidth / 2.,dy1 + posPrec[1] -YWidth / 2., dz1 + posPrec[2] - ZWidth / 2.);
                    glVertex3d(dx1 + pos[0] - XWidth / 2., dy1 + pos[1] -YWidth / 2., dz1 +pos[2] - ZWidth / 2.);
                    glEnd();
                }
            }
        }
    }

    // I define the map here to share it with the surfaces (next block)
    std::map<int, double> maxIntensityPerPoint;
    if(showHeatmap && (currentHeatmap.size() > 0)){
        if((IDheatmap < 0) || (IDheatmap >= currentHeatmap.size() + 1)){
            cerr << "ERR: Out of bounds IDheatmap" << endl;
        } else {

            for(size_t ID = 0; ID < currentHeatmap.size(); ++ID){
                if((ID == IDheatmap) || (IDheatmap == currentHeatmap.size())){
                    if(currentHeatmap[ID]){
                        size_t nH = currentHeatmap[ID]->size();
                        for(size_t i = 0; i < nH; ++i){
                            int thisPosition = (*currentHeatmap[ID])[i].first;
                            double intensity = (*currentHeatmap[ID])[i].second;
                            if(maxIntensityPerPoint.find(thisPosition) != maxIntensityPerPoint.end()){
                                intensity = max(maxIntensityPerPoint[thisPosition], intensity);
                                maxIntensityPerPoint[thisPosition] = intensity;
                            } else {
                                maxIntensityPerPoint[thisPosition] = intensity;
                            }
                        }
                    }
                }
            }
            for(std::map<int, double>::iterator it = maxIntensityPerPoint.begin(); it != maxIntensityPerPoint.end(); ++it){
                vector<int> pos = lattice::positionFromID(it->first);
                double intensity = it->second;
                glColor4d(heatmapZeroColor[0] + intensity * heatmapOneColor[0], heatmapZeroColor[1] + intensity * heatmapOneColor[1], heatmapZeroColor[2] + intensity * heatmapOneColor[2], 0.2);

                //if(intensity > 0.999) glColor4d(0.2, 1.0, 0.2, 0.2);
                if(intensity > 0.999) glColor4d(heatmapCoreColor[0], heatmapCoreColor[1], heatmapCoreColor[2], 1.0);
                glTranslated(pos[0] - XWidth / 2., pos[1] -YWidth / 2., pos[2] - ZWidth / 2.);
                glutSolidSphere(1.05*sizeSphere, 50, 50);
                glTranslated( -pos[0] + XWidth / 2.,  -pos[1] +YWidth / 2.,  -pos[2] + ZWidth / 2.);

            }
        }
    }

    int cptTouchCore = 0;
    int cptTouchHotspot = 0;
    if(visualizeInterface){
        std::map<int ,double> DistribDistancesAnyHeatmap;
        std::map<int ,double> DistribDistances100pctHeatmap;
        for(unsigned int is = 0; is < currentPoints.size(); ++is){
            vector<double> posC = currentPoints[is];

            double distanceSquareAnyHeatmap = 1000000;
            double distanceSquare100pctHeatmap = 1000000;
            // finds the closest hotspot point
            for(std::map<int, double>::iterator it = maxIntensityPerPoint.begin(); it != maxIntensityPerPoint.end(); ++it){
                vector<int> posH = lattice::positionFromID(it->first);
                vector<double> posHD = {static_cast<double>(posH[0]) - XWidth / 2., static_cast<double>(posH[1])- YWidth / 2., static_cast<double>(posH[2]) - ZWidth / 2.};
                double intensity = it->second;
                //cout << printVector(posC) << "\t" << printVector(posHD) << "\t" << sqrt(norm2(posC, posHD)) << "\n";
                distanceSquareAnyHeatmap = min(distanceSquareAnyHeatmap,sqrt(norm2(posC, posHD)));
                if(intensity > 0.999) distanceSquare100pctHeatmap = min(distanceSquare100pctHeatmap, sqrt(norm2(posC, posHD)));
            }
            if(DistribDistancesAnyHeatmap.find(floor(5*distanceSquareAnyHeatmap)) == DistribDistancesAnyHeatmap.end()){
                DistribDistancesAnyHeatmap[floor(5*distanceSquareAnyHeatmap)] = 0;
            }
            if(DistribDistances100pctHeatmap.find(floor(5*distanceSquare100pctHeatmap)) == DistribDistancesAnyHeatmap.end()){
                DistribDistances100pctHeatmap[floor(5*distanceSquare100pctHeatmap)] = 0;
            }
            DistribDistancesAnyHeatmap[floor(5*distanceSquareAnyHeatmap)] += 1;
            DistribDistances100pctHeatmap[floor(5*distanceSquare100pctHeatmap)] += 1;

            if(distanceSquareAnyHeatmap <= 1.5) {
                glColor4d(0.3, 0.5, 0.2, 1.0);
                if(distanceSquare100pctHeatmap <= 1.5) {
                    cptTouchCore++;
                    glColor4d(0.8, 0.8, 0.3, 1.0);
                } else {
                    cptTouchHotspot++;
                }
                glTranslated(posC[0], posC[1], posC[2]);
                glutSolidSphere(0.3, 50, 50);
                glTranslated(-posC[0], -posC[1], -posC[2]);
            }
        }
        cout << "For the displayed heatmaps, the minimal distances to any heatmap are" << endl;
        for(std::map<int, double>::iterator it = DistribDistancesAnyHeatmap.begin(); it != DistribDistancesAnyHeatmap.end(); ++it){
            cout << it->first/5. << "\t" << it->second << "\t[" << it->first/5. << "," << (it->first + 1) / 5. << "[\n";
        }
        cout << "For the displayed heatmaps, the minimal distances to heatmap CORES are" << endl;
        for(std::map<int, double>::iterator it = DistribDistances100pctHeatmap.begin(); it != DistribDistances100pctHeatmap.end(); ++it){
            cout << it->first/5. << "\t" << it->second << "\t[" << it->first/5. << "," << (it->first + 1) /5. << "[\n";
        }
        cout << endl;
    }


    #define sizeSphereSurface 0.12

    if(showForbPositions){
        for(unsigned int is = 0; is < currentSurfaces.size(); ++is){
            set<int>* srf = currentSurfaces[is];
            if(srf == nullptr) continue;

            //stringstream tt;
            //tt << "S="<< s->sequence;
            //drawBitmapText(tt.str().c_str(),-XWidth / 3.0, -YWidth/2.0 -is*5 , 0);

            std::set<int>::iterator it;
            for(it = srf->begin(); it != srf->end(); ++it){
                vector<int> pos = lattice::positionFromID(*it);
                //GLfloat d1[] = { 0.1, 0.4, 0.9, 1.0 };
                //glMaterialfv(GL_FRONT,GL_DIFFUSE,d1);
                //glEnable(GL_COLOR_MATERIAL);
                //glColor3f(0.1, 0.4, 0.9);
                glColor4d(0.5, 1.0, 0.3, 1.0);
                glTranslated(pos[0] - XWidth / 2., pos[1] -YWidth / 2., pos[2] - ZWidth / 2.);
                glutSolidSphere(sizeSphereSurface, 50, 50);
                glTranslated(-pos[0] + XWidth / 2., -pos[1] +YWidth / 2., -pos[2] + ZWidth / 2.);
            }
        }


        glColor4d(0.6, 0.9, 0.4, 1.0);
        glTranslated(manualPointX - XWidth / 2., manualPointY -YWidth / 2., manualPointZ - ZWidth / 2.);
        glutSolidSphere(sizeSphereSurface, 50, 50);
        glTranslated(-manualPointX + XWidth / 2., -manualPointY +YWidth / 2., -manualPointZ + ZWidth / 2.);

        for(unsigned int is = 0; is < currentPoints.size(); ++is){
            vector<double> pos = currentPoints[is];
            if(pos.size() != 3) {
                cerr << "Incorrect size of point position (currentPoints[" << is << "]\n";
                continue;
            }
            //stringstream tt;
            //tt << "S="<< s->sequence;
            //drawBitmapText(tt.str().c_str(),-XWidth / 3.0, -YWidth/2.0 -is*5 , 0);

            // Possibility to use NANs to not have lines
            if(!isnan(pos[0])){
                glColor4d(0.5, 1.0, 0.3, 1.0);
                glTranslated(pos[0] , pos[1] , pos[2] );
                glutSolidSphere(sizeSphereSurface, 50, 50);
                glTranslated(-pos[0] , -pos[1] , -pos[2] );

                if((is > 0) && (!isnan(currentPoints[is-1][0]))){
                    glBegin(GL_LINES);
                    glVertex3d(currentPoints[is-1][0], currentPoints[is-1][1], currentPoints[is-1][2]);
                    glVertex3d(pos[0], pos[1], pos[2]);
                    glEnd();
                }
            }

        }
    }

    glTranslated(-moveX, -moveY, -moveZ);

    //This is to say display everything now
    //glFlush();        // in the single buffer mode
    glutSwapBuffers();  // in the double buffers mode (see init)

    // To export the image as a bmp
    if(saving){
        int w = glutGet( GLUT_WINDOW_WIDTH );
        int h = glutGet( GLUT_WINDOW_HEIGHT );
        vector< unsigned char > buf( w * h * 3 );
        glPixelStorei( GL_PACK_ALIGNMENT, 1 );
        glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &buf[0] );

        static int cpt = 1000;
        cpt ++;
        stringstream newFile;
        newFile << PictureFileNamePrefix << "H=" << IDheatmap << "S=" << indexStructDisplayed << "im" <<  cpt << ".bmp";
        //int err =
        SOIL_save_image(newFile.str().c_str(), SOIL_SAVE_TYPE_BMP, w, h, 3, &buf[0]);

        if(visualizeInterface){
            ofstream f_wr("info_interaction.txt", ios::app);
            f_wr << PictureFileNamePrefix << "\t" << IDheatmap << "\t" << cptTouchCore << "\t" << cptTouchHotspot << "\n";
            f_wr.close();
        }
        saving = false;
    }
}


void savePicture(string fileName){
    display();
    glutPostRedisplay();
    glutSwapBuffers();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    vector< unsigned char > buf( w * h * 3 );
    glPixelStorei( GL_PACK_ALIGNMENT, 1 );

    // trick to flip the picture https://community.khronos.org/t/glreadpixels-any-possibility-to-flip-image-by-opengl/16528
    // copy starts at lower left corner of screen
//    glRasterPos2f(0,0);

//    // flip the y direction
//    glPixelZoom(1,-1);

//    // copy in place
//    glCopyPixels( 0, 0, w, h, GL_COLOR );


    glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &buf[0] );
//    vector< unsigned char > buf2( w * h * 3 );
//    for(int x = 1; x < w; ++x)
//    {
//     for (int y = 1; y < h; ++y)
//     {
//      buf2[x * w + y] = buf[x * w + (h-y)];
//     }
//    }
    //int err =
    SOIL_save_image(fileName.c_str(), SOIL_SAVE_TYPE_BMP, w, h, 3, &buf[0]);
}


/*int bytesToUsePerPixel = 24 ;  // 16 bit per channel
 *
ILuint currentImage = ilGenImage();
//ILuint ilTexName[2]; //DevIL images
//ilGenImages(2, ilTexName); //Generate DevIL image objects
//ilBindImage(ilTexName[0]);
ilBindImage(currentImage);


int sizeOfByte = sizeof( unsigned char ) ;
//int theSize = w * h * sizeOfByte * bytesToUsePerPixel ;
ilTexImage(w,h,1,bytesToUsePerPixel,GL_LUMINANCE,IL_UNSIGNED_BYTE,&buf[0]);
ilLoadL()

static int cpt = 0;
cpt ++;
QString z = QString("im") + QString::number(cpt);
ilEnable(IL_FILE_OVERWRITE);
ilSave(IL_PNG, z.toStdString().c_str());

//ILubyte *brickData, *globeData; //DevIL image data
//brickData = ilGetData(); //Get image data from DevIL

GLuint textureID;
glGenTextures(1, &textureID);
glBindTexture(GL_TEXTURE_2D, textureID);
glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE); // automatic mipmap
glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 512, 512, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
glBindTexture(GL_TEXTURE_2D, 0);*/

// create a renderbuffer object to store depth info
/*GLuint rboId;
glGenRenderbuffersEXT(1, &amp;rboId);
glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rboId);
glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, 512, 512);
glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);

// create a framebuffer object
GLuint fboId;
glGenFramebuffersEXT(1, &amp;fboId);
glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fboId);

// attach the texture to FBO color attachment point
glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, tex, 0);

// attach the renderbuffer to depth attachment point
glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, rboId);

//this is where I draw all of my scene with the exception of the sphere. too much code right here to post

// unbind FBO
glDeleteRenderbuffersEXT(1, &amp;rboId);
glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}
*/






/*
GLfloat body[][2] = { {0, 3}, {1, 1}, {5, 1}, {8, 4}, {10, 4}, {11, 5},
  {11, 11.5}, {13, 12}, {13, 13}, {10, 13.5}, {13, 14}, {13, 15}, {11, 16},
  {8, 16}, {7, 15}, {7, 13}, {8, 12}, {7, 11}, {6, 6}, {4, 3}, {3, 2},
  {1, 2} };
GLfloat arm[][2] = { {8, 10}, {9, 9}, {10, 9}, {13, 8}, {14, 9}, {16, 9},
  {15, 9.5}, {16, 10}, {15, 10}, {15.5, 11}, {14.5, 10}, {14, 11}, {14, 10},
  {13, 9}, {11, 11}, {9, 11} };
GLfloat leg[][2] = { {8, 6}, {8, 4}, {9, 3}, {9, 2}, {8, 1}, {8, 0.5}, {9, 0},
  {12, 0}, {10, 1}, {10, 2}, {12, 4}, {11, 6}, {10, 7}, {9, 7} };
GLfloat eye[][2] = { {8.75, 15}, {9, 14.7}, {9.6, 14.7}, {10.1, 15},
  {9.6, 15.25}, {9, 15.25} };
GLfloat lightZeroPosition[] = {10.0, 4.0, 10.0, 1.0};
GLfloat lightZeroColor[] = {0.8, 1.0, 0.8, 1.0}; // green-tinted
GLfloat lightOnePosition[] = {-1.0, -2.0, 1.0, 0.0};
GLfloat lightOneColor[] = {0.6, 0.3, 0.2, 1.0}; // red-tinted
GLfloat skinColor[] = {0.1, 1.0, 0.1, 1.0}, eyeColor[] = {1.0, 0.2, 0.2, 1.0};



void
extrudeSolidFromPolygon(GLfloat data[][2], unsigned int dataSize,
  GLdouble thickness, GLuint side, GLuint edge, GLuint whole)
{
  static GLUtriangulatorObj *tobj = NULL;
  GLdouble vertex[3], dx, dy, len;
  int i;
  int count = (int) (dataSize / (2 * sizeof(GLfloat)));

  if (tobj == NULL) {
    tobj = gluNewTess();  // create and initialize a GLU polygon tesselation object
    gluTessCallback(tobj, GLU_BEGIN, glBegin);
    gluTessCallback(tobj, GLU_VERTEX, glVertex2fv);  // semi-tricky
    gluTessCallback(tobj, GLU_END, glEnd);

  }
  glNewList(side, GL_COMPILE);
  glShadeModel(GL_SMOOTH);  // smooth minimizes seeing tessellation
  gluBeginPolygon(tobj);
  for (i = 0; i < count; i++) {
    vertex[0] = data[i][0];
    vertex[1] = data[i][1];
    vertex[2] = 0;
    gluTessVertex(tobj, vertex, data[i]);
  }
  gluEndPolygon(tobj);
  glEndList();
  glNewList(edge, GL_COMPILE);
  glShadeModel(GL_FLAT);  // flat shade keeps angular hands from being * * "smoothed"
  glBegin(GL_QUAD_STRIP);
  for (i = 0; i <= count; i++) {
    // mod function handles closing the edge
    glVertex3f(data[i % count][0], data[i % count][1], 0.0);
    glVertex3f(data[i % count][0], data[i % count][1], thickness);
    // Calculate a unit normal by dividing by Euclidean
       distance. We * could be lazy and use
       glEnable(GL_NORMALIZE) so we could pass in * arbitrary
       normals for a very slight performance hit.
    dx = data[(i + 1) % count][1] - data[i % count][1];
    dy = data[i % count][0] - data[(i + 1) % count][0];
    len = sqrt(dx * dx + dy * dy);
    glNormal3f(dx / len, dy / len, 0.0);
  }
  glEnd();
  glEndList();
  glNewList(whole, GL_COMPILE);
  glFrontFace(GL_CW);
  glCallList(edge);
  glNormal3f(0.0, 0.0, -1.0);  // constant normal for side
  glCallList(side);
  glPushMatrix();
  glTranslatef(0.0, 0.0, thickness);
  glFrontFace(GL_CCW);
  glNormal3f(0.0, 0.0, 1.0);  // opposite normal for other side
  glCallList(side);
  glPopMatrix();
  glEndList();
}*/


/*
void
makeDinosaur(void)
{
  extrudeSolidFromPolygon(body, sizeof(body), bodyWidth,
    BODY_SIDE, BODY_EDGE, BODY_WHOLE);
  extrudeSolidFromPolygon(arm, sizeof(arm), bodyWidth / 4,
    ARM_SIDE, ARM_EDGE, ARM_WHOLE);
  extrudeSolidFromPolygon(leg, sizeof(leg), bodyWidth / 2,
    LEG_SIDE, LEG_EDGE, LEG_WHOLE);
  extrudeSolidFromPolygon(eye, sizeof(eye), bodyWidth + 0.2,
    EYE_SIDE, EYE_EDGE, EYE_WHOLE);
  glNewList(DINOSAUR, GL_COMPILE);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor);
  glCallList(BODY_WHOLE);

  glPushMatrix();
  glTranslatef(0.0, 0.0, bodyWidth);
  glCallList(ARM_WHOLE);
  glCallList(LEG_WHOLE);
  glTranslatef(0.0, 0.0, -bodyWidth - bodyWidth / 4);
  glCallList(ARM_WHOLE);
  glTranslatef(0.0, 0.0, -bodyWidth / 4);
  glCallList(LEG_WHOLE);
  glTranslatef(0.0, 0.0, bodyWidth / 2 - 0.1);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, eyeColor);
  glCallList(EYE_WHOLE);
  glPopMatrix();
  glEndList();
}*/

#endif

// Previous ways to include openGL
/*#ifdef WINDOWS
#include "freeglut/include/GL/glut.h"
#else
#include <GL/glut.h>
#endif*/
// the openGL.dll can be found in C:/windows/SysWOW64
//file:///C:/Windows/SysWOW64/opengl32.dll
// the library can be found in:
// C:\Qt\Qt5.7.0\Tools\mingw530_32\i686-w64-mingw32\lib
// the headers are in :
// C:\Qt\Qt5.7.0\Tools\mingw530_32\i686-w64-mingw32\include
// to write text
//void showMessage(GLfloat x, GLfloat y, GLfloat z, char *message)
//{
//    glPushMatrix();
//    glDisable(GL_LIGHTING);
//    glTranslatef(x, y, z);
//    glScalef(.02, .02, .02);
//    while (*message) {
//        glutStrokeCharacter(GLUT_STROKE_ROMAN, *message);
//        message++;
//    }
//    glEnable(GL_LIGHTING);
//    glPopMatrix();
//}
// // old function, done inside display now
//void redraw(void){
//  if (newModel)
//    recalcModelView();
//  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//  //glCallList(DINOSAUR);
//  showMessage(2, 7.1, 4.1, "Spin me.");
//  glutSwapBuffers();
//}*/
