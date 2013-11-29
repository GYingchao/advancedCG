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

void pointLight(in int i, in vec3 normal, in vec3 eye, in vec3 ecPosition3)
{
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
    
   float nDotVP;       // normal . light direction
   float nDotHV;       // normal . light half vector
   float pf;           // power factor
   float attenuation;  // computed attenuation factor
   float d;            // distance from surface to light source
   vec3  VP;           // direction from surface to light position
   vec3  halfVector;   // direction of maximum highlights

   // Compute vector from surface to light position
   VP = vec3 (gl_LightSource[i].position) - ecPosition3;

   // Compute distance between surface and light position
   d = length(VP);

   // Normalize the vector from surface to light position
   VP = normalize(VP);

   // Compute attenuation
   attenuation = 1.0 ;/// (gl_LightSource[i].constantAttenuation +
      // gl_LightSource[i].linearAttenuation * d +
      // gl_LightSource[i].quadraticAttenuation * d * d);

   halfVector = normalize(VP + eye);

   nDotVP = max(0.0, dot(normal, VP));
   nDotHV = max(0.0, dot(normal, halfVector));

   pf = (nDotVP == 0.0) ? 0.0 : pow(nDotHV, gl_FrontMaterial.shininess);

   Ambient  += gl_LightSource[i].ambient * attenuation;
   Diffuse  += gl_LightSource[i].diffuse * nDotVP * attenuation;
   Specular += gl_LightSource[i].specular * pf * attenuation;
}

void main()
{
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    gl_TexCoord[0] = gl_MultiTexCoord0;    

    vec3 normal = normalize(gl_NormalMatrix * gl_Normal);

    // Eye-coordinate position of vertex, needed in various calculations
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
        
    vec4 color = gl_FrontLightModelProduct.sceneColor +
		Ambient * gl_FrontMaterial.ambient +
		Diffuse * gl_FrontMaterial.diffuse +
		Specular * gl_FrontMaterial.specular;
    color = clamp( color, 0.0, 1.0 );
    gl_FrontColor = color;
}

// Fragment shader: put it below this tag
[frag]

#version 150 compatibility

uniform sampler2D colorMap;

void main()
{   
    gl_FragColor = gl_Color * texture2D(colorMap, gl_TexCoord[0].st);
}
