// Tampere University
// TIE-52306 Computer Graphics Coding Assignment 2019
//
// Write your name and student id here:
//   Jinhyeok Yoo, 291234
//
// Mark here with an X what functionalities you implemented
// Note that different functionalities are worth different amount of points.
//
// Name of the functionality      |Done| Notes
//-------------------------------------------------------------------------------
// example functionality          | X  | Example note: control this with var YYYY
// Mandatory functionalities ----------------------------------------------------
//   Perspective projection       | X  | worked by perspective function
//   Phong shading                | X  | worked by phong function
//   Camera movement and rotation | X  | controlled by time and mouse
//   Sharp shadows                | X  | implemented by sharpshadow function but softshadow is used instead
// Extra functionalities --------------------------------------------------------
//   Attend visiting lecture 1    | X  | 
//   Attend visiting lecture 2    | X  | 
//   Tone mapping                 | X  | worked by tone_mapping function
//   PBR shading                  |    | 
//   Soft shadows                 | X  | worked by softshadow function
//   Sharp reflections            |    | 
//   Glossy reflections           |    | 
//   Refractions                  |    | 
//   Caustics                     |    | 
//   SDF Ambient Occlusions       |    | 
//   Texturing                    |    | 
//   Simple game                  |    | 
//   Progressive path tracing     |    | 
//   Basic post-processing        |    | 
//   Advanced post-processing     |    | 
//   Screen space reflections     |    | 
//   Screen space AO              |    | 
//   Simple own SDF               | X  | torus
//   Advanced own SDF             |    | 
//   Animated SDF                 | X  | twisting box
//   Other?                       |    | 


#ifdef GL_ES
precision mediump float;
#endif

#define PI 3.14159265359
#define EPSILON 0.00001

// These definitions are tweakable.

/* Minimum distance a ray must travel. Raising this value yields some performance
 * benefits for secondary rays at the cost of weird artefacts around object
 * edges.
 */
#define MIN_DIST 0.08
/* Maximum distance a ray can travel. Changing it has little to no performance
 * benefit for indoor scenes, but useful when there is nothing for the ray
 * to intersect with (such as the sky in outdoors scenes).
 */
#define MAX_DIST 20.0
/* Maximum number of steps the ray can march. High values make the image more
 * correct around object edges at the cost of performance, lower values cause
 * weird black hole-ish bending artefacts but is faster.
 */
#define MARCH_MAX_STEPS 128
/* Typically, this doesn't have to be changed. Lower values cause worse
 * performance, but make the tracing stabler around slightly incorrect distance
 * functions.
 * The current value merely helps with rounding errors.
 */
#define STEP_RATIO 0.999
/* Determines what distance is considered close enough to count as an
 * intersection. Lower values are more correct but require more steps to reach
 * the surface
 */
#define HIT_RATIO 0.001

// Resolution of the screen
uniform vec2 u_resolution;

// Mouse coordinates
uniform vec2 u_mouse;

// Time since startup, in seconds
uniform float u_time;

struct material
{
    // The color of the surface
    vec4 color;
    // You can add your own material features here!
};

// Good resource for finding more building blocks for distance functions:
// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

/* Basic box distance field.
 *
 * Parameters:
 *  p   Point for which to evaluate the distance field
 *  b   "Radius" of the box
 *
 * Returns:
 *  Distance to the box from point p.
 */
float box(vec3 p, vec3 b)
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

/* Rotates point around origin along the X axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_x(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        p.x,
        c*p.y-s*p.z,
        s*p.y+c*p.z
    );
}

/* Rotates point around origin along the Y axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_y(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        c*p.x+s*p.z,
        p.y,
        -s*p.x+c*p.z
    );
}

/* Rotates point around origin along the Z axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_z(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        c*p.x-s*p.y,
        s*p.x+c*p.y,
        p.z
    );
}

/* Each object has a distance function and a material function. The distance
 * function evaluates the distance field of the object at a given point, and
 * the material function determines the surface material at a point.
 */
float blob_distance(vec3 p)
{
    vec3 q = p - vec3(-0.5, -2.2 + abs(sin(u_time*3.0))/1.5, 2.0);
    return length(q) - 0.8 + sin(10.0*q.x)*sin(10.0*q.y)*sin(10.0*q.z)*0.07;
}

material blob_material(vec3 p)
{
    material mat;
    mat.color = vec4(1.0, 0.5, 0.3, 0.0);
    return mat;
}

float sphere_distance(vec3 p)
{
    return length(p - vec3(1.5, -1.8, 4.0)) - 1.2;
}

material sphere_material(vec3 p)
{
    material mat;
    mat.color = vec4(0.1, 0.2, 0.0, 1.0);
    return mat;
}

float room_distance(vec3 p)
{
    return max(
        -box(p-vec3(0.0,3.0,3.0), vec3(0.5, 0.5, 0.5)),
        -box(p-vec3(0.0,0.0,0.0), vec3(3.0, 3.0, 6.0))
    );
}

material room_material(vec3 p)
{
    material mat;
    mat.color = vec4(1.0, 1.0, 1.0, 1.0);
    if(p.x <= -2.98) mat.color.rgb = vec3(1.0, 0.0, 0.0);
    else if(p.x >= 2.98) mat.color.rgb = vec3(0.0, 1.0, 0.0);
    return mat;
}

float crate_distance(vec3 p)
{
    return box(rot_y(p-vec3(-1,-1,5), u_time), vec3(1, 2, 1));
}

material crate_material(vec3 p)
{
    material mat;
    mat.color = vec4(1.0, 1.0, 1.0, 1.0);

    vec3 q = rot_y(p-vec3(-1,-1,5), u_time) * 0.98;
    if(fract(q.x + floor(q.y*2.0) * 0.5 + floor(q.z*2.0) * 0.5) < 0.5)
    {
        mat.color.rgb = vec3(0.0, 1.0, 1.0);
    }
    return mat;
}

float torus_distance(vec3 p)
{
    p.x -= 1.5;
    p.y -= 1.0;
    vec2 q = vec2(length(p.xy) - 0.7, p.z - 6.0);
    return length(q) - 0.3;
}

material torus_material(vec3 p)
{
    material mat;
    mat.color.rgb = vec3(0.534,0.318,1.000);
    return mat;
}

float twistBox_distance(vec3 p)
{
    p -= vec3(1.5, -2.4, 1.4);
    vec3 q = rot_y(p, sin(u_time)*1.5*p.y);
    return box(q , vec3(0.6));
}

material twistBox_material(vec3 p)
{
    material mat;
    mat.color.rgb = vec3(0.720,0.720,0.000);
    return mat;
}

/* The distance function collecting all others.
 *
 * Parameters:
 *  p   The point for which to find the nearest surface
 *  mat The material of the nearest surface
 *
 * Returns:
 *  The distance to the nearest surface.
 */
float map(
    in vec3 p,
    out material mat
){
    float min_dist = MAX_DIST*2.0;
    float dist = 0.0;

    dist = blob_distance(p);
    if(dist < min_dist) {
        mat = blob_material(p);
        min_dist = dist;
    }

    dist = room_distance(p);
    if(dist < min_dist) {
        mat = room_material(p);
        min_dist = dist;
    }

    dist = crate_distance(p);
    if(dist < min_dist) {
        mat = crate_material(p);
        min_dist = dist;
    }

    dist = sphere_distance(p);
    if(dist < min_dist) {
        mat = sphere_material(p);
        min_dist = dist;
    }

    // Add your own objects here!
    dist = torus_distance(p);
    if(dist < min_dist) {
        mat = torus_material(p);
        min_dist = dist;
    }

    dist = twistBox_distance(p);
    if(dist < min_dist) {
        mat = twistBox_material(p);
        min_dist = dist;
    }

    return min_dist;
}

/* Calculates the normal of the surface closest to point p.
 *
 * Parameters:
 *  p   The point where the normal should be calculated
 *  mat The material information, produced as a byproduct
 *
 * Returns:
 *  The normal of the surface.
 *
 * See http://www.iquilezles.org/www/articles/normalsSDF/normalsSDF.htm if
 * you're interested in how this works.
 */
vec3 normal(vec3 p, out material mat)
{
    const vec2 k = vec2(1.0, -1.0);
    return normalize(
        k.xyy * map(p + k.xyy * EPSILON, mat) +
        k.yyx * map(p + k.yyx * EPSILON, mat) +
        k.yxy * map(p + k.yxy * EPSILON, mat) +
        k.xxx * map(p + k.xxx * EPSILON, mat)
    );
}

/* Finds the closest intersection of the ray with the scene.
 *
 * Parameters:
 *  o           Origin of the ray
 *  v           Direction of the ray
 *  max_dist    Maximum distance the ray can travel. Usually MAX_DIST.
 *  p           Location of the intersection
 *  n           Normal of the surface at the intersection point
 *  mat         Material of the intersected surface
 *  inside      Whether we are marching inside an object or not. Useful for
 *              refractions.
 *
 * Returns:
 *  true if a surface was hit, false otherwise.
 */
bool intersect(
    in vec3 o,
    in vec3 v,
    in float max_dist,
    out vec3 p,
    out vec3 n,
    out material mat,
    bool inside
) {
    float t = MIN_DIST;
    float dir = inside ? -1.0 : 1.0;
    bool hit = false;

    for (int i = 0; i < MARCH_MAX_STEPS; ++i)
    {
        p = o + t * v;
        float dist = dir * map(p, mat);
        
        hit = abs(dist) < HIT_RATIO * t;

        if(hit || t > max_dist) break;

        t += dist * STEP_RATIO;
    }

    n = normal(p, mat);

    return hit;
}

/* Finds the point where the shadow should be cast.
 *
 * Parameters:
 *  p           Point of the surface
 *  ld          Direction of the light
 *  max_dist    Maximum distance the ray can travel
 *  mat         Material of the intersected surface
 *
 * Returns:
 *  0.25 if the shadow should be cast, 1.0 otherwise.
 */
float sharpshadow(vec3 p, vec3 ld, float max_dist, material mat)
{
    float t = MIN_DIST;
    
    for(int i = 0; i < MARCH_MAX_STEPS; i++)
    {
		float dist = map(p + t * ld, mat);
        
        if (abs(dist) < HIT_RATIO * t)	
            return 0.25;
        
        t += dist * STEP_RATIO;
        
        if (t > max_dist)
            break;
    }
    return 1.0;
}

/* Finds the point where the shadow should be cast.
 * Reference: http://iquilezles.org/www/articles/rmshadows/rmshadows.htm
 *
 * Parameters:
 *  p           Point of the surface
 *  ld          Direction of the light
 *  max_dist    Maximum distance the ray can travel
 *  mat         Material of the intersected surface
 *
 * Returns:
 *  Intensity of the shadow.
 */
float softshadow(vec3 p, vec3 ld, float max_dist, material mat)
{
	float res = 1.0;
    float t = MIN_DIST;
    
    for(int i = 0; i < MARCH_MAX_STEPS; i++)
    {
		float dist = map(p + t * ld, mat);
        res = min(res, 5.0*dist/t);
        t += dist * STEP_RATIO;
        
        if (res < EPSILON || t > max_dist)
            break;   
    }
    return clamp(res, 0.25, 1.0);
}

/* Creates the perspective projection
 * from exercise 3
 *
 * Parameters:
 *  ro   origin of the view ray
 *  rt   destination of the view ray
 *  q    screen coordinates
 *
 * Returns:
 *  Direction of the view ray.
 */
vec3 perspective(vec3 ro, vec3 rt, vec2 q){
    // We make a camera coordinate system.
    vec3 z = normalize(rt-ro);			// Direction vector from camera to the vanishing point
    vec3 y = vec3(0.0, 1.0, 0.0);		// Up direction vector
    vec3 x = normalize(cross(y, z));	// left direction vector perpendicular to y and z
    // The perspective cameara can be implemeted
    // by setting the direction of the ray as screen coordinates plus z-value weights 'radians(80.0)'.
    // The closer to the center (0,0) of the screen coordinates,
    // the closer the direction of ray is to the perfect forward direction vector (0,0,1).
    // So we implement the perspective camera using this way in the ray tracer.
    // Plus, by multiplying the camera coordinate system matrix (by transforming the coordinate system),
    // we create perspective regardless of the position and direction of the view.
    // And if the origin of the view ray changes, the camera is rotated to look at the vanishing point.
    return normalize(mat3(x, y, z)*vec3(q, radians(80.0)));
}

/* Calculates the phong reflection model
 * from exercise 3
 *
 * Parameters:
 *  n    surface normal
 *  v    view direction
 *  l    light direction
 *  mat  Material of the objects
 *
 * Returns:
 *  Color of the pixel.
 */
vec3 phong(vec3 n, vec3 v, vec3 l, in material mat)
{
    // Reflection coefficient
	const float Ka = 0.3;
	const float Kd = 1.0;
	const float Ks = 1.0;
	const float ShininessVal = 60.0;
    
    // Lambert's cosine law
    float lambertian = max(dot(n, l), 0.0);
    float specular = 0.0;
    if(lambertian > 0.0) {
        vec3 r = reflect(-l, n);      // Reflected light vector
        // Compute the specular term
        float specAngle = max(dot(r, v), 0.0);
        specular = pow(specAngle, ShininessVal);
    }
    return Ka * mat.color.rgb + Kd * lambertian * mat.color.rgb + Ks * specular * vec3(1.000,1.000,1.000);
}

/* Calculates the color of the pixel, based on view ray origin and direction.
 *
 * Parameters:
 *  o   Origin of the view ray
 *  v   Direction of the view ray
 *
 * Returns:
 *  Color of the pixel.
 */
vec3 render(vec3 o, vec3 v)
{   
    // This lamp is positioned at the hole in the roof.
	vec3 lampPos = vec3(0.0, 3.1, 3.0);
    
    vec3 p, n;
    material mat;

    // Compute intersection point along the view ray.
    intersect(o, v, MAX_DIST, p, n, mat, false);

    // Add some lighting code here!
    mat.color.rgb = phong(n, normalize(o - p), normalize(lampPos - p), mat)
        		 // * sharpshadow(p, normalize(lampPos - p), 0.8, mat);
        			* softshadow(p, normalize(lampPos - p), 0.8, mat);

    return mat.color.rgb;
}

/* Adjusts the gamma and exposure for the color of the pixel.
 *
 * Parameters:
 *  color   pixel color to apply tone mapping
 *
 * Returns:
 *  Tone-mapped pixel color.
 */
vec3 tone_mapping(vec3 color)
{
    float gamma = 1.5;
	float exposure = 1.0;
	color = clamp(exposure * color, 0.0, 1.0);
	color = pow(color, vec3(1.0 / gamma));
	return color;
}

void main()
{
    // This is the position of the pixel in normalized device coordinates.
    vec2 uv = (gl_FragCoord.xy/u_resolution)*2.0-1.0;
    // Calculate aspect ratio
    float aspect = u_resolution.x/u_resolution.y;
    // The position of mouse cursor
    vec2 mouse = (u_mouse/u_resolution)*2.0-1.0;

    vec3 ro = vec3(-mouse, 3.0*abs(sin(u_time/2.0))-5.0);
    vec3 rt = vec3(0.0, 0.0, 0.0);

    vec3 rd = perspective(ro, rt, vec2(uv.x * aspect, uv.y));
    
    vec3 color = render(ro, rd);
    gl_FragColor = vec4(tone_mapping(color), 1.0);
}
