#ifdef WIN32
	#include <windows.h>
#endif
#include "GL/glew.h"
#include "GL/glut.h"

#include "vector3.h"
#include "model_obj.h"
#include "bitmap.h"
#include "glShader.h"

#include <map>


// Enumeration
enum EnumDisplayMode { WIREFRAME = 0, HIDDENLINE, FLATSHADED, SMOOTHSHADED, TEXTURESMOOTHSHADED, SHADERGOURAUD, SHADERPHONG,/*SHADOWMAP, SHADOWVOLUME,*/ MODENUM };

char* g_DisplayModeNames[] = {
	"Wireframe",
	"Hidden Line", 
	"Flat Shaded", 
	"Smooth Shaded", 
	"Textured Smooth Shaded", 
	"Custom Shader - Gouraud Shading (Per Vert)", 
	"Custom Shader - Phong Shading (Per Frag)", 
//	"w/ Shadow Map", 
//	"w/ Shadow Volume"
};

typedef std::map<std::string, GLuint> ModelTextures;

// variables
EnumDisplayMode displayMode = TEXTURESMOOTHSHADED;	// current display mode
int mainMenu, displayMenu;		// glut menu handlers
int winWidth, winHeight;		// window width and height
double winAspect;				// winWidth / winHeight;
int lastX, lastY;				// last mouse motion position
bool leftDown, middleDown, middleUp, shiftDown;		// mouse down and shift down flags
float sphi = 90.0, stheta = 45.0, sdepth = 10;	// for simple trackball
float xpan = 0.0, ypan = 0.0;				// for simple trackball
float zNear = 1.0, zFar = 100.0;
float g_fov = 45.0;
Vector3f g_center;
float g_fFrameTime = 0;
int g_iFrameCount = 0;
float g_fFPS = 0;
ModelOBJ g_model;	// OBJ mesh representation
ModelTextures       g_modelTextures;
GLuint		g_nullTexture = 0;
GLShader	g_shaderPerVertLight;
GLShader	g_shaderPerFragLight;

float				g_maxAnisotrophy = 1.0f;
bool                g_enableTextures = true;

// functions
void SetBoundingBox();
void InitGL();
void InitMenu();
void InitGeometry();
void MenuCallback(int value);
void ReshapeFunc(int width, int height);
void DisplayFunc();
void IdleFunc();
void DrawWireframe();
void DrawHiddenLine();
void DrawFlatShaded();
void DrawSmoothShaded();
void DrawTexturedSmoothShaded();
void DrawShaderPerVertexLighting();
void DrawShaderPerFragmentLighting();

void KeyboardFunc(unsigned char ch, int x, int y);
void MouseFunc(int button, int state, int x, int y);
void MotionFunc(int x, int y);
void CalculateFPS();
void DrawText(float x, float y, char *string);

GLuint LoadTexture(const char *pszFilename);
GLuint CreateNullTexture(int width, int height);
void LoadModel(const char *pszFilename);
void UnloadModel();

void SetBoundingBox() {
	
	double PI = 3.14159265358979323846;
	float radius = g_model.getRadius();
	g_model.getCenter(g_center.X(), g_center.Y(), g_center.Z());
	zNear    = float(0.2 * radius / sin(0.5 * g_fov * PI / 180.0) );
	zFar     = zNear + 2.0f * radius;
	sdepth = zNear + radius;
	zNear *= 0.1f;
	zFar *= 10.f;
}

// init openGL environment
void InitGL() {
	GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f}; 
	GLfloat light0Position[] = { 3.0f, 10.0f, 0.0f, 1.0f }; 
	GLfloat light0Ambient[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	GLfloat light0Diffuse[] = { 1.0f, 1.0f, 0.0f, 1.0f };
	GLfloat light0Specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

	GLfloat light1Position[] = { 3.0f, -10.0f, 0.0f, 1.0f }; 
	GLfloat light1Ambient[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	GLfloat light1Diffuse[] = { 1.0f, 0.0f, 0.0f, 1.0f };
	GLfloat light1Specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(800, 600);
	glutCreateWindow("COMP541 Mesh Viewer");

	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		throw std::runtime_error("GLEW initialization error.\n");
		exit(0);
	}

	fprintf( stdout, "OpenGL Version = '%s'\n", glGetString(GL_VERSION) );
	fprintf( stdout, "Shading Language Version = '%s'\n", glGetString(GL_SHADING_LANGUAGE_VERSION) );

	if (!glewIsSupported("GL_VERSION_3_2"))
	{
		fprintf(stderr, "Warning: Your graphics card or the driver may be old.");
		fprintf(stderr, "    The next assignment will require OpenGL 3.2 and GLSL 1.50 support.");
		fprintf(stderr, "    Please change the \"#version\" tag to 110 in your GLSL shader file.");
	}
	else if (!glewIsSupported("GL_VERSION_2_0"))
	{
		fprintf(stderr, "Warning: Your graphics card or the driver may be too old.");
		fprintf(stderr, "    This assignment requires at least OpenGL 2.0 and GLSL 1.10 support.");
		throw std::runtime_error("GLEW initialization error.\n");
	}

	glEnable(GL_MULTISAMPLE);
	glEnable(GL_ARB_compatibility);
	glClearColor(0.3f, 0.5f, 0.9f, 1.0f);
	glPolygonOffset(1.0, 1.0);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glDisable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_LIGHTING);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);
	glLightfv (GL_LIGHT0, GL_POSITION, light0Position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0Ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0Diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0Specular);
	glEnable(GL_LIGHT0);

	glLightfv (GL_LIGHT1, GL_POSITION, light1Position);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light1Ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1Diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light1Specular);
	glEnable(GL_LIGHT1);

	// Set texture filtering options
	glActiveTexture(GL_TEXTURE1);
	glEnable(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	// load per-vertex lighting shader
	const char* filenamevert = "..\\shaders\\blinn_phong_vert.glsl";

	fprintf(stdout, "Compiling shader \"%s\".\n", filenamevert);
	std::string infoLog;
	g_shaderPerVertLight.LoadShaderProgramFromFile(filenamevert, infoLog);
	if(!g_shaderPerVertLight.GetShader())
	{
		fprintf(stderr, "%s\n", infoLog.c_str());
		throw std::runtime_error("Failed to load Gouraud shader.\n" + infoLog);
	}
	else
		fprintf(stdout, "Shader (Gouraud: per-vertex shading) compiled successfully.\n");

	// load per-fragment lighting shader
	const char* filenamefrag = "..\\shaders\\blinn_phong_frag.glsl";

	fprintf(stdout, "Compiling shader \"%s\".\n", filenamefrag);
	g_shaderPerFragLight.LoadShaderProgramFromFile(filenamefrag, infoLog);
	if(!g_shaderPerFragLight.GetShader())
	{
		fprintf(stderr, "%s\n", infoLog.c_str());
		throw std::runtime_error("Failed to load Phong shader.\n" + infoLog);
	}
	else
		fprintf(stdout, "Shader (Phong: per-fragment shading) compiled successfully.\n");

	// Create null texture for consistently shading models without a texture
	g_nullTexture = CreateNullTexture(2, 2);

	// Set callback functions
	glutReshapeFunc(ReshapeFunc);
	glutDisplayFunc(DisplayFunc);
	glutIdleFunc(IdleFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
}

// init right-click menu
void InitMenu() {
	displayMenu = glutCreateMenu(MenuCallback);
	for(int i = 0; i < MODENUM; ++i)
	{
		glutAddMenuEntry(g_DisplayModeNames[i], i);
	}

	mainMenu = glutCreateMenu(MenuCallback);
	glutAddSubMenu("Display", displayMenu);
	glutAddMenuEntry("Exit", 99);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void ChangeDisplayMode(EnumDisplayMode mode)
{
	if(int(mode) >= MODENUM || int(mode) < 0 )
		return;
	displayMode = mode;
	fprintf(stdout, "Switching to %s mode.\n", g_DisplayModeNames[int(mode)]);
	glutPostRedisplay();
}

// GLUT menu callback function
void MenuCallback(int value) {
	switch (value) {
	case 99: exit(0); break;
	default: 
		ChangeDisplayMode(EnumDisplayMode(value));
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
	// Set the projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(g_fov, winAspect, zNear, zFar);

	// Set the model*view matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
	glTranslatef(xpan, ypan, -sdepth);
	glRotatef(-stheta, 1.0, 0.0, 0.0);
	glRotatef(sphi, 0.0, 0.0, 1.0);
	glTranslatef(-g_center[0], -g_center[1], -g_center[2]);

	// clear the framebuffer and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	switch (displayMode) {
	case WIREFRAME: DrawWireframe(); break;
	case HIDDENLINE: DrawHiddenLine(); break;
	case FLATSHADED: DrawFlatShaded(); break;
	case SMOOTHSHADED: DrawSmoothShaded(); break;
	case TEXTURESMOOTHSHADED: DrawTexturedSmoothShaded(); break;
	case SHADERGOURAUD: DrawShaderPerVertexLighting(); break;
	case SHADERPHONG: DrawShaderPerFragmentLighting(); break;
	}

	//  Print the FPS to the window
	char strBuf[100];
	sprintf_s(strBuf, 100, "FPS: %4.1f", g_fFPS);
	DrawText(-0.9f, -0.9f, strBuf);

	glutSwapBuffers();
}


void DrawModelShaded()
{
	const ModelOBJ::Mesh *pMesh = 0;
	const ModelOBJ::Material *pMaterial = 0;
	const ModelOBJ::Vertex *pVertices = 0;
	ModelTextures::const_iterator iter;
	GLuint texture = 0;
	
	// Iterate all the object meshes in the OBJ file
	for (int i = 0; i < g_model.getNumberOfMeshes(); ++i)
	{
		pMesh = &g_model.getMesh(i);
		pMaterial = pMesh->pMaterial;
		pVertices = g_model.getVertexBuffer();

		// Set mesh-specific material properties
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, pMaterial->ambient);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, pMaterial->diffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, pMaterial->specular);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, pMaterial->shininess * 128.0f);

		// Default texture is all white
		texture = g_nullTexture;

		// Apply the correct texture if exists
		if (g_enableTextures)
		{
			iter = g_modelTextures.find(pMaterial->colorMapFilename);

			if (iter != g_modelTextures.end())
				texture = iter->second;
		}

		glActiveTexture(GL_TEXTURE0);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texture);

		// Bind vertex position buffer
		if (g_model.hasPositions())
		{
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, g_model.getVertexSize(),
				g_model.getVertexBuffer()->position);
		}

		// Bind texture coordinate buffer
		if (g_model.hasTextureCoords())
		{
			glClientActiveTexture(GL_TEXTURE0);
			glEnableClientState(GL_TEXTURE_COORD_ARRAY);
			glTexCoordPointer(2, GL_FLOAT, g_model.getVertexSize(),
				g_model.getVertexBuffer()->texCoord);
		}

		// Bind normal direction buffer
		if (g_model.hasNormals())
		{
			glEnableClientState(GL_NORMAL_ARRAY);
			glNormalPointer(GL_FLOAT, g_model.getVertexSize(),
				g_model.getVertexBuffer()->normal);
		}

		// Bind tangent direction buffer
		if (g_model.hasTangents())
		{
			glClientActiveTexture(GL_TEXTURE1);
			glEnableClientState(GL_TEXTURE_COORD_ARRAY);
			glTexCoordPointer(4, GL_FLOAT, g_model.getVertexSize(),
				g_model.getVertexBuffer()->tangent);
		}

		// Draw all the triangles in one batch. Yay!
		glDrawElements(GL_TRIANGLES, pMesh->triangleCount * 3, GL_UNSIGNED_INT,
			g_model.getIndexBuffer() + pMesh->startIndex);

		// Unbind the input buffers
		if (g_model.hasTangents())
		{
			glClientActiveTexture(GL_TEXTURE1);
			glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		}
		if (g_model.hasNormals())
			glDisableClientState(GL_NORMAL_ARRAY);
		if (g_model.hasTextureCoords())
		{
			glClientActiveTexture(GL_TEXTURE0);
			glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		}
		if (g_model.hasPositions())
			glDisableClientState(GL_VERTEX_ARRAY);
	}

	glBindTexture(GL_TEXTURE_2D, 0);
}

// Wireframe render function
void DrawWireframe() {
	g_enableTextures = false;
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDisable(GL_LIGHTING);
	glColor3f(0, 0, 0);
	glDisable(GL_CULL_FACE);
	DrawModelShaded();
	glEnable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

// Hidden Line render function
void DrawHiddenLine() {
	g_enableTextures = false;
	glShadeModel(GL_FLAT); 
	glEnable(GL_POLYGON_OFFSET_FILL);
	glDrawBuffer(GL_NONE);	// depth only pass, prime the depth buffer.
	DrawModelShaded();
	glDrawBuffer(GL_BACK);
	glDisable(GL_LIGHTING);
	glColor3f(0, 0, 0);
	glDisable(GL_POLYGON_OFFSET_FILL);
	DrawWireframe();
	glEnable(GL_LIGHTING);
}

// Flat Shaded render function
void DrawFlatShaded() {
	g_enableTextures = false;
	glShadeModel(GL_FLAT); 
	glEnable(GL_LIGHTING);
	DrawModelShaded();
}


// Smooth Shaded render function
void DrawSmoothShaded() { 
	g_enableTextures = false;
	glShadeModel(GL_SMOOTH); 
	glEnable(GL_LIGHTING);
	DrawModelShaded();
}

// Textured Smooth Shaded render function
void DrawTexturedSmoothShaded() { 
	g_enableTextures = true;
	glShadeModel(GL_SMOOTH); 
	glEnable(GL_LIGHTING);
	DrawModelShaded();
}

// Render using the programmable shader for per-vertex lighting
void DrawShaderPerVertexLighting() { 
	g_enableTextures = true;
	glUseProgram(g_shaderPerVertLight.GetShader());

	// Update shader parameters.
	glUniform1i(glGetUniformLocation(
		g_shaderPerVertLight.GetShader(), "colorMap"), 0);

	glUniform1f(glGetUniformLocation(
		g_shaderPerVertLight.GetShader(), "g_fFrameTime"), g_fFrameTime);

	DrawModelShaded();

	glUseProgram(0);
}

// Render using the programmable shader for per-fragment lighting
void DrawShaderPerFragmentLighting() { 
	g_enableTextures = true;
	glUseProgram(g_shaderPerFragLight.GetShader());

	// Update shader parameters.
	glUniform1i(glGetUniformLocation(
		g_shaderPerFragLight.GetShader(), "colorMap"), 0);

	glUniform1f(glGetUniformLocation(
		g_shaderPerFragLight.GetShader(), "g_fFrameTime"), g_fFrameTime);

	DrawModelShaded();

	glUseProgram(0);
}


// GLUT keyboard callback function
void KeyboardFunc(unsigned char ch, int x, int y) {
	switch (ch) {
	case '1': case '2': case '3': case '4': case '5': 
	case '6': case '7': case '8': case '9':  
		ChangeDisplayMode(EnumDisplayMode(ch - '1'));
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
}

// GLUT mouse motion callback function
void MotionFunc(int x, int y) {
	if (leftDown)
		if(!shiftDown) { // rotate
			sphi += (float)(x - lastX) / 4.0f;
			stheta += (float)(lastY - y) / 4.0f;
		} else { // pan
			xpan += (float)(x - lastX)*sdepth/zNear/winWidth;
			ypan += (float)(lastY - y)*sdepth/zNear/winHeight;
		}
	// scale
	if (middleDown) sdepth += (float)(lastY - y) / 50.0f;

	lastX = x;
	lastY = y;
	glutPostRedisplay();
}

void IdleFunc() {
	// Get the number of milliseconds since glutInit was called
	g_fFrameTime = (float)glutGet(GLUT_ELAPSED_TIME); 

	CalculateFPS();
	//  Call display function (draw the current frame)
	glutPostRedisplay ();
}

void CalculateFPS()
{
	static float previousTime;
	//  Increase frame count
	g_iFrameCount++;

	//  Calculate time passed
	int timeInterval = int(g_fFrameTime - previousTime);

	// Do it in every second
	if(timeInterval > 1000)
	{
		//  calculate the number of frames per second
		g_fFPS = g_iFrameCount / (timeInterval / 1000.0f);
		previousTime = g_fFrameTime;
		g_iFrameCount = 0;
	}
}

// main function
void main(int argc, char **argv) {
	glutInit(&argc, argv);
	InitGL();
	InitMenu();
	if (argc>=2) 
		LoadModel(argv[1]);
	else
	{
		fprintf(stderr, "Error: No input model specified.\n");
		fprintf(stderr, "Usage: pa3.exe ..\\models\\venus.obj.\n");
		throw std::runtime_error("No input model specified.\n");
		exit(0);
	}
	SetBoundingBox();

	glutMainLoop();
}



GLuint LoadTexture(const char *pszFilename)
{
	GLuint id = 0;
	Bitmap bitmap;

	if (bitmap.loadPicture(pszFilename))
	{
		// The Bitmap class loads images and orients them top-down.
		// OpenGL expects bitmap images to be oriented bottom-up.
		bitmap.flipVertical();

		glGenTextures(1, &id);
		glBindTexture(GL_TEXTURE_2D, id);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

		if (g_maxAnisotrophy > 1.0f)
		{
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT,
				g_maxAnisotrophy);
		}

		gluBuild2DMipmaps(GL_TEXTURE_2D, 4, bitmap.width, bitmap.height,
			GL_BGRA_EXT, GL_UNSIGNED_BYTE, bitmap.getPixels());
	}

	return id;
}

GLuint CreateNullTexture(int width, int height)
{
	// Create an empty white texture. This texture is applied to OBJ models
	// that don't have any texture maps. This trick allows the same shader to
	// be used to draw the OBJ model with and without textures applied.

	int pitch = ((width * 32 + 31) & ~31) >> 3; // align to 4-byte boundaries
	std::vector<GLubyte> pixels(pitch * height, 255);
	GLuint texture = 0;

	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_BGRA,
		GL_UNSIGNED_BYTE, &pixels[0]);

	return texture;
}

void LoadModel(const char *pszFilename)
{
	// Import the OBJ file and normalize to unit length.

	SetCursor(LoadCursor(0, IDC_WAIT));
	
	fprintf(stdout, "Loading model \"%s\". \n", pszFilename);

	if (!g_model.import(pszFilename))
	{
		SetCursor(LoadCursor(0, IDC_ARROW));
		throw std::runtime_error("Failed to load model.");
		exit(0);
	}

	g_model.normalize();

	// Load any associated textures.
	// Note the path where the textures are assumed to be located.

	const ModelOBJ::Material *pMaterial = 0;
	GLuint textureId = 0;
	std::string::size_type offset = 0;
	std::string filename;

	fprintf(stdout, "Loading materials. \n");

	for (int i = 0; i < g_model.getNumberOfMaterials(); ++i)
	{
		pMaterial = &g_model.getMaterial(i);

		// Look for and load any diffuse color map textures.

		if (pMaterial->colorMapFilename.empty())
			continue;

		// Try load the texture using the path in the .MTL file.
		textureId = LoadTexture(pMaterial->colorMapFilename.c_str());

		if (!textureId)
		{
			offset = pMaterial->colorMapFilename.find_last_of('\\');

			if (offset != std::string::npos)
				filename = pMaterial->colorMapFilename.substr(++offset);
			else
				filename = pMaterial->colorMapFilename;

			// Try loading the texture from the same directory as the OBJ file.
			textureId = LoadTexture((g_model.getPath() + filename).c_str());
		}

		if (textureId)
			g_modelTextures[pMaterial->colorMapFilename] = textureId;

		// Look for and load any normal map textures.

		if (pMaterial->bumpMapFilename.empty())
			continue;

		// Try load the texture using the path in the .MTL file.
		textureId = LoadTexture(pMaterial->bumpMapFilename.c_str());

		if (!textureId)
		{
			offset = pMaterial->bumpMapFilename.find_last_of('\\');

			if (offset != std::string::npos)
				filename = pMaterial->bumpMapFilename.substr(++offset);
			else
				filename = pMaterial->bumpMapFilename;

			// Try loading the texture from the same directory as the OBJ file.
			textureId = LoadTexture((g_model.getPath() + filename).c_str());
		}

		if (textureId)
			g_modelTextures[pMaterial->bumpMapFilename] = textureId;
	}

	fprintf(stdout, "Model loading completed. \n");

	SetCursor(LoadCursor(0, IDC_ARROW));
}

void UnloadModel()
{
	SetCursor(LoadCursor(0, IDC_WAIT));

	ModelTextures::iterator i = g_modelTextures.begin();

	while (i != g_modelTextures.end())
	{
		glDeleteTextures(1, &i->second);
		++i;
	}

	g_modelTextures.clear();
	g_model.destroy();

	SetCursor(LoadCursor(0, IDC_ARROW));
}

//  Draws a string at the specified coordinates.
void DrawText(float x, float y, char *string)
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -1.0f);
	glRasterPos2f(x, y);

	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	for (int i = 0, len = (int)strlen(string); i < len; i++)
	{ 
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, (int)string[i]);
	}
	glPopAttrib();
}
