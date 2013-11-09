//Per-fragment Blinn-Phong shader (Phong shading) for two single directional light sources.

[vert]

#version 150 compatibility

// Time since the first frame (in milliseconds) 
uniform float g_fFrameTime;

// data to be passed down to a later stage
out vec3 normal;
out vec4 ecPosition;

void main()
{
    // Compute View(eye) space surface normal at the current vertex
    normal = normalize(gl_NormalMatrix * gl_Normal);
	
	// TODO(4): add an non-trivial interesting animation to the model 
	//       by changing "vAnimatedPos" using a function of "g_fFrameTime"
	vec4 vAnimatedPos = gl_Vertex;

    // Eye-coordinate position of vertex, needed in lighting computation
    ecPosition = gl_ModelViewMatrix * vAnimatedPos;

    // Compute the projection space position of the current vertex
    gl_Position = gl_ModelViewProjectionMatrix * vAnimatedPos;
    // Pass on te texture coordinate
    gl_TexCoord[0] = gl_MultiTexCoord0;    
}

[frag]

#version 150 compatibility

uniform sampler2D colorMap;
uniform float materialAlpha;

// data passed down and interpolated from the vertex shader
in vec3 normal;
in vec4 ecPosition;

// global variables used in auxilary functions
vec4 Ambient;
vec4 Diffuse;
vec4 Specular;

void pointLight(in int i, in vec3 normal, in vec3 eye, in vec3 ecPosition3)
{
    // TODO(3): Copy the completed per-vertex point light computation to here
	
	//First we compute the ambient constribution
    Ambient  += gl_LightSource[i].ambient;

	//Then we compute the diffuse lighting
	//	Assume default light source position is in eye space.
	//vec4 ls_ecPosition = gl_ModelViewMatrix * gl_LightSource[i].position;
	//vec3 ls_ecPosition3 = (vec3(ls_ecPosition)) / ls_ecPosition.w;
	vec3 ls_ecPosition3 = (vec3(gl_LightSource[i].position)) / gl_LightSource[i].position.w;
	vec3 L = normalize(ls_ecPosition3 - ecPosition3);
	float dot_product = dot(normal, L);
	if(dot_product < 0.0) dot_product = 0.0;
    Diffuse  += gl_LightSource[i].diffuse * dot_product;

	//Then we compute the specular lighting component
	float specDot;
	//	Visibility test (Whether the angle between N and L are more than 90 degree)
	if (dot(normal, L) < 0.0) specDot = 0.0;
	else {
		vec3 V = normalize(eye - ecPosition3);
		vec3 H = normalize(L + V);
		specDot = dot(normal, H);
		if(specDot <= 0.0) specDot = 0.0;
		else specDot = pow(specDot, gl_FrontMaterial.shininess);
	} 
    Specular += gl_LightSource[i].specular * specDot;
}


void main()
{   
    // Process passed-down interpolated attributes
    vec3 n = normalize(normal);
    vec3 ecPosition3 = (vec3 (ecPosition)) / ecPosition.w;
    
    //TODO(3): Move the Blinn-Phong computation from the vertex shader to the pixel shader.
    //      Important difference: Perform the texture lookup and modulate by the 
    //      lighting result *prior to* adding the specular component.

	// Clear the light intensity accumulators
    Ambient  = vec4 (0.0);
    Diffuse  = vec4 (0.0);
    Specular = vec4 (0.0);
    vec3 eye = vec3 (0.0, 0.0, 1.0);
    
    // Compute point light contributions
    pointLight(0, normal, eye, ecPosition3);
    pointLight(1, normal, eye, ecPosition3);
        
    // TODO(1): Add ambient, diffuse and specular contributions to equation below.
   vec4 color = gl_FrontLightModelProduct.sceneColor + gl_FrontMaterial.ambient*Ambient + gl_FrontMaterial.diffuse*Diffuse + gl_FrontMaterial.specular*Specular;
        
	// Clamp color to [0, 1]
    color = clamp( color, 0.0, 1.0 );
    gl_FragColor = color;
}
