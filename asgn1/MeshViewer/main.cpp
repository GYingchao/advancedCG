#ifdef WIN32
	#include <windows.h>
#endif
#include <GL/glut.h>
#include "mesh.h"
#include "matrix.h"
#include "OpenGLProjector.h"
extern Weight currentWeight;
enum EnumDisplayMode { WIREFRAME, HIDDENLINE, FLATSHADED, SMOOTHSHADED, COLORSMOOTHSHADED };

enum Mode 
{ 
	Viewing, 
	Selection 
};
Mode currentMode = Viewing;

// variables
int displayMode = FLATSHADED;	// current display mode
int mainMenu, displayMenu;		// glut menu handlers
int winWidth, winHeight;		// window width and height
double winAspect;				// winWidth / winHeight;
int lastX, lastY;				// last mouse motion position
int currSelectedVertex = -1;         // current selected vertex
bool leftDown, middleDown, middleUp, shiftDown;		// mouse down and shift down flags
double sphi = 90.0, stheta = 45.0, sdepth = 10;	// for simple trackball
double xpan = 0.0, ypan = 0.0;				// for simple trackball
double zNear = 1.0, zFar = 100.0;
double g_fov = 45.0;
Vector3d g_center;
double g_sdepth;
Mesh mesh;	// our mesh

// functions
void SetBoundaryBox(const Vector3d & bmin, const Vector3d & bmax);
void InitGL();
void InitMenu();
void InitGeometry();
void MenuCallback(int value);
void ReshapeFunc(int width, int height);
void DisplayFunc();
void DrawWireframe();
void DrawHiddenLine();
void DrawFlatShaded();
void DrawSmoothShaded();
void DrawColorSmoothShaded();
void DrawSelectedVertices();

void KeyboardFunc(unsigned char ch, int x, int y);
void MouseFunc(int button, int state, int x, int y);
void MotionFunc(int x, int y);
void SelectVertexByPoint();
void DeleteSelectedVertex(int vertex);

// help funcs
void DrawBoundaryEdges();

void SetBoundaryBox(const Vector3d & bmin, const Vector3d & bmax) {
	//double PI = 3.14159265358979323846;
	double radius = bmax.Distance(bmin);
	g_center = 0.5 * (bmin+bmax);
	zNear    = 0.2 * radius / sin(0.5 * g_fov * PI / 180.0);
	zFar     = zNear + 2.0 * radius;
	g_sdepth = zNear + radius;
	zNear *= 0.1;
	zFar *= 10;
	sdepth = g_sdepth;
}

// init openGL environment
void InitGL() {
	GLfloat light0Position[] = { 0, 1, 0, 1.0 }; 

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
	glutCreateWindow("Comp541 Mesh Viewer");
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glPolygonOffset(1.0, 1.0);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glLightfv (GL_LIGHT0, GL_POSITION, light0Position);
	glEnable(GL_LIGHT0);

	glutReshapeFunc(ReshapeFunc);
	glutDisplayFunc(DisplayFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
}

// init right-click menu
void InitMenu() {
	displayMenu = glutCreateMenu(MenuCallback);
	glutAddMenuEntry("Wireframe", WIREFRAME);
	glutAddMenuEntry("Hidden Line", HIDDENLINE);
	glutAddMenuEntry("Flat Shaded", FLATSHADED);
	glutAddMenuEntry("Smooth Shaded", SMOOTHSHADED);
	glutAddMenuEntry("Color Smooth Shaded", COLORSMOOTHSHADED);

	mainMenu = glutCreateMenu(MenuCallback);
	glutAddSubMenu("Display", displayMenu);
	glutAddMenuEntry("Exit", 99);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

// init geometry (if no input argument is provided)
void InitGeometry() {
	const int VSIZE = 4;
	const int HESIZE = 12;
	const int FSIZE = 4;
	int i;
	Vertex *v[VSIZE];
	HEdge *he[HESIZE];
	Face *f[FSIZE];
	
	for (i=0; i<VSIZE; i++) {
		v[i] = new Vertex();
		mesh.vList.push_back(v[i]);
	}
	v[0]->SetPosition(Vector3d(0.0, 0.0, 0.0));
	v[1]->SetPosition(Vector3d(10.0, 0.0, 0.0));
	v[2]->SetPosition(Vector3d(0.0, 10.0, 0.0));
	v[3]->SetPosition(Vector3d(0.0, 0.0, 10.0));

	v[0]->SetNormal(Vector3d(-0.577, -0.577, -0.577));
	v[1]->SetNormal(Vector3d(0.0, -0.7, -0.7));
	v[2]->SetNormal(Vector3d(-0.7, 0.0, -0.7));
	v[3]->SetNormal(Vector3d(-0.7, -0.7, 0.0));

	for (i=0; i<FSIZE; i++) {
		f[i] = new Face();
		mesh.fList.push_back(f[i]);
	}

	for (i=0; i<HESIZE; i++) {
		he[i] = new HEdge();
		mesh.heList.push_back(he[i]);
	}
	for (i=0; i<FSIZE; i++) {
		int base = i*3;
		SetPrevNext(he[base], he[base+1]);
		SetPrevNext(he[base+1], he[base+2]);
		SetPrevNext(he[base+2], he[base]);
		SetFace(f[i], he[base]);
	}
	SetTwin(he[0], he[4]);
	SetTwin(he[1], he[7]);
	SetTwin(he[2], he[10]);
	SetTwin(he[3], he[8]);
	SetTwin(he[5], he[9]);
	SetTwin(he[6], he[11]);
	he[0]->SetStart(v[1]); he[1]->SetStart(v[2]); he[2]->SetStart(v[3]);
	he[3]->SetStart(v[0]); he[4]->SetStart(v[2]); he[5]->SetStart(v[1]);
	he[6]->SetStart(v[0]); he[7]->SetStart(v[3]); he[8]->SetStart(v[2]);
	he[9]->SetStart(v[0]); he[10]->SetStart(v[1]); he[11]->SetStart(v[3]);
	v[0]->SetHalfEdge(he[3]);
	v[1]->SetHalfEdge(he[0]);
	v[2]->SetHalfEdge(he[1]);
	v[3]->SetHalfEdge(he[2]);
}

// GLUT menu callback function
void MenuCallback(int value) {
	switch (value) {
	case 99: exit(0); break;
	default: 
		displayMode = value;
		glutPostRedisplay();
		break;
	}
}

// GLUT reshape callback function
void ReshapeFunc(int width, int height) {
	winWidth = width;
	winHeight = height;
	winAspect = (double)width/(double)height;
	glViewport(0, 0, width, height);
	glutPostRedisplay();
}

// GLUT display callback function
void DisplayFunc() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(g_fov, winAspect, zNear, zFar);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
	glTranslatef(xpan, ypan, -sdepth);
	glRotatef(-stheta, 1.0, 0.0, 0.0);
	glRotatef(sphi, 0.0, 1.0, 0.0);
	glTranslatef(-g_center[0], -g_center[1], -g_center[2]);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	switch (displayMode) {
	case WIREFRAME: DrawWireframe(); break;
	case HIDDENLINE: DrawHiddenLine(); break;
	case FLATSHADED: DrawFlatShaded(); break;
	case SMOOTHSHADED: DrawSmoothShaded(); break;
	case COLORSMOOTHSHADED: DrawColorSmoothShaded(); break;
	}
	//mesh.ComputeVertexNormals();
	DrawSelectedVertices();
	DrawBoundaryEdges();
	
	glutSwapBuffers();
}

// Wireframe render function
void DrawWireframe() {
	HEdgeList heList = mesh.Edges();
	HEdgeList bheList = mesh.BoundaryEdges();
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINES);
	size_t i;
	for (i=0; i<heList.size(); i++) {
		glVertex3dv(heList[i]->Start()->Position().ToArray());
		glVertex3dv(heList[i]->End()->Position().ToArray());
	}
	glColor3f(1, 0, 0);
	for (i=0; i<bheList.size(); i++) {
		glVertex3dv(bheList[i]->Start()->Position().ToArray());
		glVertex3dv(bheList[i]->End()->Position().ToArray());
	}
	glEnd();
}

// Hidden Line render function
void DrawHiddenLine() {
	FaceList fList = mesh.Faces();
	glShadeModel(GL_FLAT); 
	glEnable(GL_POLYGON_OFFSET_FILL);
	glColor3f(0, 0, 0);
	glBegin(GL_TRIANGLES);
	for (size_t i=0; i<fList.size(); i++) {
		Face *f = fList[i];
		const Vector3d & pos1 = f->HalfEdge()->Start()->Position();
		const Vector3d & pos2 = f->HalfEdge()->End()->Position();
		const Vector3d & pos3 = f->HalfEdge()->Next()->End()->Position();
		glVertex3dv(pos1.ToArray());
		glVertex3dv(pos2.ToArray());
		glVertex3dv(pos3.ToArray());
	}
	glEnd();
	glDisable(GL_POLYGON_OFFSET_FILL);

	DrawWireframe();
}

// Flat Shaded render function
void DrawFlatShaded() {
	FaceList fList = mesh.Faces();
	glShadeModel(GL_FLAT); 
	glEnable(GL_LIGHTING);
	glColor3f(0.4f, 0.4f, 1.0f);
	glBegin(GL_TRIANGLES) ;
	for (size_t i=0; i<fList.size(); i++) {
		Face *f = fList[i];
		const Vector3d & pos1 = f->HalfEdge()->Start()->Position();
		const Vector3d & pos2 = f->HalfEdge()->End()->Position();
		const Vector3d & pos3 = f->HalfEdge()->Next()->End()->Position();
		Vector3d normal = (pos2-pos1).Cross(pos3-pos1);
		normal /= normal.L2Norm();
		glNormal3dv(normal.ToArray());
		glVertex3dv(pos1.ToArray());
		glVertex3dv(pos2.ToArray());
		glVertex3dv(pos3.ToArray());
	}
	glEnd();
	glDisable(GL_LIGHTING);
}

// Smooth Shaded render function
void DrawSmoothShaded() { 
	FaceList fList = mesh.Faces();
	glShadeModel(GL_SMOOTH); 
	glEnable(GL_LIGHTING);
	glColor3f(0.4f, 0.4f, 1.0f);
	glBegin(GL_TRIANGLES) ;
	for (size_t i=0; i<fList.size(); i++) {
		Face *f = fList[i];
		Vertex * v1 = f->HalfEdge()->Start();
		Vertex * v2 = f->HalfEdge()->End();
		Vertex * v3 = f->HalfEdge()->Next()->End();
		glNormal3dv(v1->Normal().ToArray());
		glVertex3dv(v1->Position().ToArray());
		glNormal3dv(v2->Normal().ToArray());
		glVertex3dv(v2->Position().ToArray());
		glNormal3dv(v3->Normal().ToArray());
		glVertex3dv(v3->Position().ToArray());
	}
	glEnd();
	glDisable(GL_LIGHTING);
}

void DrawColorSmoothShaded() {
	FaceList fList = mesh.Faces();
	glShadeModel(GL_SMOOTH); 
	glEnable(GL_LIGHTING);
	glColor3f(0.4f, 0.4f, 1.0f);
	glBegin(GL_TRIANGLES) ;
	for (size_t i=0; i<fList.size(); i++) {
		Face *f = fList[i];
		Vertex * v1 = f->HalfEdge()->Start();
		Vertex * v2 = f->HalfEdge()->End();
		Vertex * v3 = f->HalfEdge()->Next()->End();
		glNormal3dv(v1->Normal().ToArray());
		glColor3dv(v1->Color().ToArray());
		//glColor3f(1.0f, 1.0f, 0.0f);
		glVertex3dv(v1->Position().ToArray());

		glNormal3dv(v2->Normal().ToArray());
		glColor3dv(v2->Color().ToArray());
		//glColor3f(1.0f, 1.0f, 0.0f);
		glVertex3dv(v2->Position().ToArray());

		glNormal3dv(v3->Normal().ToArray());
		glColor3dv(v3->Color().ToArray());
		//glColor3f(1.0f, 1.0f, 0.0f);
		glVertex3dv(v3->Position().ToArray());
	}
	glEnd();
	glDisable(GL_LIGHTING);
}


// draw the selected ROI vertices on the mesh
void DrawSelectedVertices()
{
	VertexList vList = mesh.Vertices();
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	size_t i;
	for (i=0; i<vList.size(); i++) 
	{
		if (vList[i]->Flag())
		{
			switch (vList[i]->Flag() % 3)
			{
			case 0: // handle vertices
				glColor3f(1.0, 0.3, 0.3);
				break;
			case 1: 
				glColor3f(1.0, 0.0, 0.0); 
				break;
			case 2: 
				glColor3f(0.3, 0.3, 1.0); 
				break;
			}
			glVertex3dv(vList[i]->Position().ToArray());
		}
	}
	glEnd();
}

// GLUT keyboard callback function
void KeyboardFunc(unsigned char ch, int x, int y) {
	switch (ch) {
	case 'c':
	case 'C':
		if(currentWeight == Uniform) {
			currentWeight = Cotangent;
		} else {
			currentWeight = Uniform;
		}
		break;
	case 'u':
	case 'U':
		/************************************************************************/
		/* activate the following code if you finish the corresponding functions*/
		cout << "Do the UmbrellaSmoothing.." << endl;
 		mesh.UmbrellaSmooth();
 		mesh.ComputeVertexNormals();
 		mesh.ComputeVertexCurvatures();
		/************************************************************************/
		break;
	case 's':
	case 'S':
		/************************************************************************/
		/* activate the following code if you finish the corresponding functions*/
		cout << "Do the Implicit Umbrella Smoothing.." << endl;
 		mesh.ImplicitUmbrellaSmooth();
 		mesh.ComputeVertexNormals();
 		mesh.ComputeVertexCurvatures();
		/************************************************************************/
		break;
	case '1':	// key '1'
		currentMode = Viewing;
		cout << "Viewing mode" << endl;
		break;
	case '2':	// key '2'
		currentMode = Selection;
		cout << "Selection mode" << endl;
		break;
	case 'D':
	case 'd':
		if( currentMode == Selection){
			cout << "vertex deleted" << endl;
			VertexList vList = mesh.Vertices();
			Vertex* v = vList[currSelectedVertex];

			// for non-boundary selection
			if(!v->IsBoundary()){
				OneRingHEdge ring (v);
				HEdge *cure = NULL;
				HEdge* B1 = v->HalfEdge()->Twin()->Prev();// previously created new boundary edge
				HEdge* B2 = NULL;// old edge at the currtent new boundary edge postion
				HEdge* B3 = NULL;// firstly created boundary edge 
				unsigned int temp_ct = 0;
				while (cure=ring.NextHEdge()) { 
					B2 = cure->Next();
					// only to change IsBoundary
					HEdge* newB = NULL;
					newB = new HEdge(true);
					newB->SetFlag(B2->Flag());
					newB->SetStart(B2->Start());
					newB->SetValid(B2->IsValid());
						
					SetPrevNext(B2->Prev(),newB);
					SetPrevNext(newB, B2->Next());
					SetTwin(B2->Twin(), newB);
					
					newB->Start()->SetHalfEdge(newB);
					mesh.bheList.push_back(newB); 
					for(size_t i =0; i<mesh.heList.size(); i++){ // erase the old half edges 
						if( mesh.heList[i] == B2){//mesh.heList[i] == cure || mesh.heList[i] == cure->Prev() ||
							mesh.heList.erase(mesh.heList.begin()+i);
							delete B2;
							B2 = NULL;
							break;
						}
					}
				}
				OneRingHEdge ring2examine (v); // check the vertices around the target v
				cure = NULL;
				while (cure=ring2examine.NextHEdge()){
				
					// find vertices connected only to boundary edges and delete them
					Vertex* adjv = cure->Twin()->Start();
					//B1 = B1->Next();
					OneRingHEdge adjRing (adjv);
					HEdge * he = NULL;
					bool removeAdjV = true;
					HEdge* prevE=NULL;
					HEdge* nextE=NULL;
					while (he=adjRing.NextHEdge()){
						bool heB = he->IsBoundary();
						bool heTwinB = he->Twin()->IsBoundary();
						if(heB^heTwinB){ // exclusive or
							if (heB == true){ // ending
								if(!(he->Prev()->IsBoundary() && !he->Prev()->Twin()->IsBoundary())) //already connected
									nextE = he;
							}
							else { // starting
								if(!(he->Twin()->Next()->IsBoundary() && !he->Twin()->Next()->Twin()->IsBoundary()))
									prevE = he->Twin();
							}
							he->Start()->SetHalfEdge(he);
							removeAdjV=false;
						}
						else if (heB == false && heTwinB == false) // triangle not affected 
						{
							he->Start()->SetHalfEdge(he);
							if(he->Twin()->Start() != v)
								removeAdjV=false;
						}
						else if (heB == true && heTwinB == true){ 
						}
					}
					if(removeAdjV){

						if(prevE != NULL || nextE != NULL){
							cout<<"error: attention 9"<<endl;
							exit(0);
						}
						cout<<"one vertex removed"<<endl;
						for(size_t i=0;i<mesh.vList.size();i++){
							
							if(mesh.vList[i] == adjv){mesh.vList.erase(mesh.vList.begin()+i); delete adjv; break;}
						}
						//delete adjv;
						//adjv= NULL;
					}
					else{
						if(prevE == NULL || nextE == NULL){
							cout<<"attention 10, this should appear only once"<<endl;
						}else{
							SetPrevNext(prevE,nextE);
						}
					}
				}

				// remove the spoke of the 'wheel'
				OneRingHEdge ring2delete (v);
				cure = NULL;
				while (cure=ring2delete.NextHEdge()){
					for(size_t i=0; i<mesh.fList.size();i++){ // erase the triangle face
						Face * tempf= cure->LeftFace();
						if(tempf!=NULL && mesh.fList[i] == tempf){ 
							mesh.fList.erase(mesh.fList.begin()+i); delete tempf;  break;}
					}
					for(size_t i =0; i<mesh.heList.size(); i++){
						if(mesh.heList[i]==cure|| mesh.heList[i]==cure->Prev()){
							mesh.heList.erase(mesh.heList.begin()+i);
							i--;
						}
					}
					for(size_t i =0; i<mesh.bheList.size(); i++){
						if(mesh.bheList[i]==cure|| mesh.bheList[i]==cure->Prev()){
							mesh.bheList.erase(mesh.bheList.begin()+i);
							i--;
						}
					}
					//delete cure->Prev();
					//delete cure;
				}
				

				for(size_t i =0; i<mesh.bheList.size(); i++){
					if(mesh.bheList[i]->IsBoundary() && mesh.bheList[i]->Twin()->IsBoundary()){
						mesh.bheList.erase(mesh.bheList.begin()+i);
						i--;
					}
				}
			}

			// for boundary selection
			else{
				OneRingHEdge ring (v);
				HEdge *cure = NULL;
				HEdge* B1 = v->HalfEdge()->Twin()->Prev();// previously created new boundary edge
				HEdge* B2 = NULL;// old edge at the currtent new boundary edge postion
				HEdge* B3 = NULL;// firstly created boundary edge 
				unsigned int temp_ct = 0;
				while (cure=ring.NextHEdge()) { 
					B2 = cure->Next();
					//B2->SetPrev(B1);
					//B2->SetNext(B2->Next()->Twin()->Next());
					
					// only to change IsBoundary
					HEdge* newB = NULL;
					bool finalLink = true;
					if(!B2->IsBoundary()&& !(B2->Twin()->IsBoundary())){ // case 1: cure far away from boundary
						
						newB = new HEdge(true);
						newB->SetFlag(B2->Flag());
						newB->SetStart(B2->Start());
						newB->SetValid(B2->IsValid());
						
						SetPrevNext(B1,newB);
						SetTwin(B2->Twin(), newB);
						
						// the original hedge could be removed later
						// so restore it here
						newB->Start()->SetHalfEdge(newB); 

						B1 = newB;
						mesh.bheList.push_back(newB); // add the new boundary edge
						//delete B2;
					}
					else if(B2->IsBoundary()){ // case 2: edge is on boundary
						if(temp_ct == 0) finalLink=false;
						newB = B2;
						SetPrevNext(B1,newB); // this B1 needs taken care of if this case happen when temp_ct ==0
						//SetTwin(B2->Twin(), newB);
						newB->Start()->SetHalfEdge(newB);
						B1 = cure->Prev()->Prev();
						if(!(B1->IsBoundary() && !B1->Twin()->IsBoundary())){
							cout<<"error: attention 002"<<endl;
						}
					}
					else{
						cout<<"error: attention 006"<<endl; // not handled
					}
					if(temp_ct == 0 ){
						if(finalLink)
							B3 = B1;
						else{
							B3 = newB; // start in case 2
						}
					} // firstly created boundary edge
					temp_ct++;
					
					for(size_t i=0; i<mesh.fList.size();i++){ // erase the triangle face
						if(cure->LeftFace()!=NULL && mesh.fList[i] == cure->LeftFace()){ 
							mesh.fList.erase(mesh.fList.begin()+i); break;}
					}
					for(size_t i =0; i<mesh.heList.size(); i++){ // erase the old half edges (a triangle)
						if(mesh.heList[i] == cure || mesh.heList[i] == cure->Prev() || mesh.heList[i] == cure->Next()){
							mesh.heList.erase(mesh.heList.begin()+i);
							i--;
						}
					}
					for(size_t i =0; i<mesh.bheList.size(); i++){ // erase the old half edges (a triangle)
						if(mesh.bheList[i] == cure || mesh.bheList[i] == cure->Prev()){
							mesh.bheList.erase(mesh.bheList.begin()+i);
							i--;
						}
					}
				}
				
				SetPrevNext(B1, B3); // connect the finaly and firstly created boundary edges
				
			}
			// erase the vertex
			cout<<"selected vertex removed"<<endl;
			for(size_t i=0;i<mesh.vList.size();i++){
				if(mesh.vList[i] == v){mesh.vList.erase(mesh.vList.begin()+i); break;}
			}
		}
		break;
	case 'm':
	case 'M':
		mesh.DisplayMeshInfo();
		break;
	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}

// GLUT mouse callback function
void MouseFunc(int button, int state, int x, int y) {
	
	lastX = x;
	lastY = y;
	leftDown = (button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN);
	middleDown = (button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN);
	middleUp = (button == GLUT_MIDDLE_BUTTON) && (state == GLUT_UP);
	shiftDown = (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
	
	if (currentMode == Selection && state == GLUT_UP)
	{
	    if (middleUp)
		{
		    if (currSelectedVertex != -1)
		    {
		        mesh.Vertices()[currSelectedVertex]->SetFlag(0);
		        currSelectedVertex = -1;
		    }
		}
		else
	    {
		    SelectVertexByPoint();
		}

		lastX = lastY = 0;
		glutPostRedisplay();
	}
}

// GLUT mouse motion callback function
void MotionFunc(int x, int y) {
	if (leftDown)
		if(!shiftDown) { // rotate
			sphi += (double)(x - lastX) / 4.0;
			stheta += (double)(lastY - y) / 4.0;
		} else { // pan
			xpan += (double)(x - lastX)*sdepth/zNear/winWidth;
			ypan += (double)(lastY - y)*sdepth/zNear/winHeight;
		}
	// scale
	if (middleDown) sdepth += (double)(lastY - y) / 10.0;

	lastX = x;
	lastY = y;
	glutPostRedisplay();
}


// select a mesh point
void SelectVertexByPoint()
{
	// get the selection position
	int x = lastX, y = winHeight - lastY;
	Vector3d u(x,y,0);
	
	OpenGLProjector projector;

	VertexList vList = mesh.Vertices();

    double mindis = 1e6; int selectedIndex = -1;
	for (size_t i=0; i<vList.size(); i++) 
	{
		Vector3d v = projector.Project(vList[i]->Position());
		
		vList[i]->SetFlag(0); // set back to unselected
		
		if (projector.GetDepthValue((int)v.X(), (int)v.Y()) - v.Z() <= -1e-4)
		{
		    continue;
		}
		v.Z() = 0;
		
		double dist = (u - v).L2Norm();
		if (dist < mindis)
		{
		    mindis = dist;
		    selectedIndex = (int)i;
		}
	}
	
	if (selectedIndex != -1)
	{
	    currSelectedVertex = selectedIndex;
	    vList[selectedIndex]->SetFlag(1);
	}
	
}

// help functions
void DrawBoundaryEdges(){

	glColor3f(1.0, 0.0, 0.0);
	glLineWidth(5.0);
	
	for(size_t i=0; i<mesh.bheList.size();i++){
		glBegin(GL_LINES);
		glVertex3dv(mesh.bheList[i]->Start()->Position().ToArray());
		glVertex3dv(mesh.bheList[i]->End()->Position().ToArray());
		glEnd();
		//cout<<*(mesh.bheList[i]->Start()->Position().ToArray())<<endl;
	}
}

// main function
void main(int argc, char **argv) { // add arg[] in Project -> pa1 properties -> configeration properties -> debugging
	glutInit(&argc, argv);
	InitGL();
	InitMenu();
	if (argc>=2) {
		cout<<"loading obj file "<<argv[1]<<endl;
		
		char str[80];
		strcpy(str, "models/");
		strcat(str, argv[1]);
		mesh.LoadObjFile(str);
	}
	else InitGeometry();
	SetBoundaryBox(mesh.MinCoord(), mesh.MaxCoord());

    mesh.DisplayMeshInfo();
	
	/************************************************************************/
	/* activate the following code if you finish the corresponding functions*/
 	mesh.ComputeVertexNormals();

	mesh.ComputeVertexCurvatures();
	/************************************************************************/

	glutMainLoop();
}


