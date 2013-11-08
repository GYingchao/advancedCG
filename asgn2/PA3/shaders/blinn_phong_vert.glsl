//Per-vertex Blinn-Phong shader (Goroud shading) for two single directional light sources.

// Vertex shader: put it below this tag
[vert]

// Shader code version: 1.5 with backward compatibility mode
// change to 110 if your graphics board does not support GLSL 1.5
#version 150 compatibility

// global variables shared between functions
vec4 Ambient;
vec4 Diffuse;
vec4 Specular;

// Built-in variables and states that are useful:
//
// Light attributes:
// gl_LightSource[i].position	: vec4	: position of lightsource i (same below)
// gl_LightSource[i].ambient	: vec4	: ambient contribution
// gl_LightSource[i].diffuse	: vec4	: diffuse contribution
// gl_LightSource[i].specular	: vec4	: specular contribution
//
// Material attributes:
// gl_FrontMaterial.shininess	: float : specular exponential term of the current shaded material (same below)
// gl_FrontMaterial.ambient     : vec4  : ambient reflective factor
// gl_FrontMaterial.diffuse     : vec4  : diffuse reflective factor
// gl_FrontMaterial.specular    : vec4  : specular reflective factor
// gl_FrontLightModelProduct.sceneColor : vec4 : the product of the ambient material property 
//                          for front-facing surfaces and the global ambient light for the scene.


void pointLight(in int i, in vec3 normal, in vec3 eye, in vec3 ecPosition3)
{
    // TODO(1): compute the ambient, diffuse and specular contributions of 
    // point light i, using the parameters and the built-in variables above
    // Parameters: 
    //   i: light index (1/2 in this assignment)
    //   normal: eye(view) space surface normal at the shaded vertex
    //   eye: eye(view) space camera position
    //   ecPosition3: eye(view) space position of the shaded vertex
    Ambient  += 0;
    Diffuse  += 0;
    Specular += 0;
}

void main()
{
    // Compute the projection space position of the current vertex
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    
    // Pass on te texture coordinate
    gl_TexCoord[0] = gl_MultiTexCoord0;    

    // Compute View(eye) space surface normal at the current vertex
    vec3 normal = normalize(gl_NormalMatrix * gl_Normal);

    // Eye-coordinate (view space) position of vertex, needed in lighting computation
    vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
    vec3 ecPosition3 = (vec3 (ecPosition)) / ecPosition.w;

	// Clear the light intensity accumulators
    Ambient  = vec4 (0.0);
    Diffuse  = vec4 (0.0);
    Specular = vec4 (0.0);
    vec3 eye = vec3 (0.0, 0.0, 1.0);
    
    // Compute point light contributions
    pointLight(0, normal, eye, ecPosition3);
    pointLight(1, normal, eye, ecPosition3);
        
    // TODO(1): Add ambient, diffuse and specular contributions to equation below.
    vec4 color = gl_FrontLightModelProduct.sceneColor + 0;
        
	// Clamp color to [0, 1]
    color = clamp( color, 0.0, 1.0 );
    
    // Output color
    gl_FrontColor = color;
}

// Fragment shader: put it below this tag
[frag]

#version 150 compatibility

// Texture sampler of the current surface
uniform sampler2D colorMap;

void main()
{   
    // Output color directly
    // TODO(2): sample the texture color from "colorMap" using the 
    //       passed down texture coordinate gl_TexCoord[0], and
    //       modulate the interpolated lighting color by it.
    gl_FragColor = gl_Color;
}