---
title: \textbf{Smoothed Particle Hydrodynamics for Real-Time Fluid Simulation for Unity}
author: "Yohib Hussain"
fontsize: 12pt
mainFont: "Times New Roman"
geometry: margin=1in
header-includes:
  - \usepackage{graphicx}
---

\hrule

## Abstract

Realistically animated fluids can add substantial realism to interactive applications such as virtual surgery simulators or computer games. In this paper we try to implement a Fluid Simulation using Smooth Particle Hydrodynamics (SPH) to simulate fluids with free surfaces in Unity Game Engine. The method is an extension of the SPH-based technique by Desbrun to animate highly deformable bodies adapting it for interactive real-time use in Unity.

My implementation leverages compute shaders to effciently parallelize particle interactions on the GPU, enabling simulations with thousands of particles at interactive frame rates. Key physical effects such as pressure, velocity and external forces are incorporated, and rendering techniques are applied to visulize fluid dynamics.

Results demonstrate that the method produces visually plausible fluid motion with good performance across different particle counts. While not physically exact compared to high-resolution grid solvers, the approach balances realism and effciency, making it suitable for applications where interactivity is critical. Future improvements may include surface tension modeling, multi-phase fluids, and integration with rigid-body dynamics.

\vspace{1em}

\hrule

## Introduction

### Motivation

Fluids (i.e liquids and gases) play an important role in every day life. Examples of fluid phenomena are wind, weather, ocean waves, waves induced by ships or simply pouring of a glass of water. As simple and ordinary these phenomena may seem, as complex and difficult it is to simulate them. Even though Computational Fluid Dynamics (CFD) is a well established research area with a long history, there are still many open research problems in the field. The reason for the complexity of fluid behaviour is the complex interplay of various phenomena such as convection, diffusion, turbulence and surface tension. Fluid phenomena are typically simulated off-time and then visualized in a second step e.g. in aerodynamics or optimization of turbines or pipes with the goal of being accurate as possible.

\newpage

Less accurate methods that allow the simulation of fluid effects in real-time open up a variety of new applications. In the fields mentioned above real-time methods help to test whether a certain concept is promising during the design phase. Other applications for real-time simulation techniques for fluids are medical simulators, computer games or any type of virtual environment. We will take a look into the computer game side of the fluid simulation

### Topic

The focus of this thesis is Smoothed Particle Hydrodynamics (SPH), a particle-based method for simulating the flow of fluids. While SPH can be also applied to gases, in this work the emphasis will be place on liquid fluids only.

## Smooth Particle Hydrodynamics

Although Smoothed Particle Hydrodynamics (SPH) was developed for the simulation of astrophysical problems, the method in general enough to be used in any kind of fluid simulation.

SPH is an interpolation method for particle systems. With SPH, fleid quantities that are only defined at the discrete particle locations can be evaluted anywhere in space. For this purpose, SPH distributes quantities in a local neighborhood of each particle using radial symmetrical smoothing kernels.

Accoring to SPH, a scalar quantity $\textit{A}$ is interpolated at location $\textbf{r}$ by a weighted sum of contributions from all particles:

\begin{align}
A_{s}(\textbf{r}) = \sum_{j}m_j\frac{A_j}{\rho_j}W(\textbf{r}-\textbf{r}_j,h)
\end{align}

where $\textit{j}$ iterates over all the particles, $m_j$ is the mass of the particle \textit{j}, $\textbf{r}_j$ its position, $\rho_j$ is the density and $\textit{A}_j$ is the field quantity at $\textbf{r}_j$.

The function $W(r,h)$ is called the smoothing kernel with core radius $h$. It defines the effective range of interaction in the simulation: We can use this to get the estimates on the density, pressure and viscosity of the given particle.

The particle mass and density appear in Eqn. $(1)$ because each particle $i$ represents a certain volume $V_i = \frac{m_i}{\rho_i}$. While the mass $m_i$ is constant throughout the simulation and, in our case, the same for all the particles, the density $\rho_i$ varies and needs to be evaluated at every time step. Through substitution into Eqn $(1)$ we get for the density at location $\textbf{r}$

\begin{align}
\rho_S(\textbf{r}) = \sum_jm_j\frac{\rho_j}{\rho_j}W(\textbf{r} - \textbf{r}_j,h)  = \sum_jm_jW(\textbf{r} - \textbf{r}_j,h)
\end{align}


## Modelling Fluids with Particles

In the Eulerian (grid based) formulation, isothermal fluids are described by a velocity field $v$, a density field $\rho$ and a pressure field $p$. The evolution of these quentities over time is given by two equations. The first equation assures convservation of mass.

equation over here.
\begin{align}
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \textbf{v}) = 0
\end{align}

while the Navier-Stokes equation formulates conservation of momentum

equation over here.
\begin{align}
\rho\left( \frac{\partial \textbf{v}}{\partial t} + \textbf{v} \cdot \nabla \textbf{v} \right) = - \nabla p + \rho \textbf{g} + \mu \nabla^{2} \textbf{v}
\end{align}

where $g$ is the external force density field and $\mu$ the viscosity of the fluid. Many forms of the Navier-Stokes equation appear in the literature. Eqn. $(4)$ represents a simplified version for incompressible fluids.

In the SPH formulation, the accleration of a particle $i$ is obtained from the Newton's second law:

\begin{align}
\frac{d \textbf{v}_i}{dt} = \frac{\textbf{F}_i}{m_i}
\end{align}

where $\textbf{v}_i$ is the velocity of the particle and $i$ and $\textbf{f}_i$ and $\rho_i$ are the force density field and the density field evaluated at the location of particle $i$, repectively.

This formulation is equivalent to the continuum expression

\begin{align}
\frac{d \textbf{v}_i}{dt} = \frac{\textbf{f}_i}{\rho_i}
\end{align}

where $\textbf{f}_i$ denotes the force per unity volume and $\rho_i$ the density, as often written in SPH literature. 

Now we will see the force density terms using SPH.

\clearpage

### Pressure

Application of the SPH rule described in Eqn $(1)$ to the pressure term $-\nabla p$ yields

\begin{align}
\textbf{f}^{pressure}_i = -\nabla p(\textbf{r}_i) = -\sum_j m_j \frac{p_j}{\rho_j} \nabla W(\textbf{r}_i - \textbf{r}_j , h)
\end{align}

Unfortunately, this force is not symmetric as can be seen when only two partciles interact. Since the gradient of the kernel is zero as its center, particle $i$ only uses the pressure of particle $j$ to compute its pressure force and vice versa. Because the pressures at the locations of the two particles are not equal in general, the pressure forces will not be symmertic. Different ways of symmertrization of Eqn. $(7)$ have been proposed in the literature. But I have went with this

\begin{align}
\textbf{f}^{pressure}_i = -\nabla p(\textbf{r}_i) = -\sum_j m_j \frac{p_i + p_j}{\rho_j} \nabla W(\textbf{r}_i - \textbf{r}_j , h)
\end{align}

The so computed pressure force is symmetric because it uses the arithmetic mean of the pressures of interacting particles.

Since particles only carry the three quantities mass, position and velocity, the pressure at particle locations has to be evaluated first. This is done in two steps. Eqn $(2)$ yields the density at the location of the particle. Then, the pressure can be computed via the ideal gas state equation

\begin{align}
p = k \rho
\end{align}

where $k$ is the gas constant that depends on the temperature. In my simulation I have used a modified version of Eqn $(9)$

\begin{align}
p = k (\rho - \rho_0)
\end{align}

where $\rho_0$ is the rest density. Which I have given as a constant

### External Forces

In the Unity Project which supports gravity, collision forces and forces caused by user interaction. These forces are applied directly to the particles without the use of SPH. When the particles collide with any solid object such as the blue box in my Unity Simulation, we simply push them out of the object and reflect the velocity component that is perpendicular to the object's surface.

\clearpage

## Smoothing kernels

Smoothing kernels are fundamental in SPH as they determine how particle properties are interpolated over space. They provide stability, accuracy and computational effciency to the simulation, such as density estimation, pressure forces, and viscosity. I have used two types of smoothing kernels

* Poly6 Kernel: This kernel is used for density estimation. It's smooth shape ensures that conributions from nearby particles are weighted appropriately, avoiding sharp discontinuities.

\begin{align}
W_{\text{poly6}}(\mathbf{r}, h) = \frac{315}{64\pi h^9} 
\begin{cases} 
(h^2 - r^2)^3 & 0 \le r \le h \\ 
0 & \text{otherwise} 
\end{cases}
\end{align}

* Spiky Kernel: Used for calculating pressure and viscosity forces. Its gradient is sharper near the particle center, which improves force calculations and help maintain stability.

\begin{align}
W_{\text{spiky}}(\mathbf{r}, h) = \frac{15}{\pi h^6} 
\begin{cases} 
(h - r)^3 & 0 \le r \le h \\ 
0 & \text{otherwise} 
\end{cases}
\end{align}

\begin{figure}[ht!]
    \centering
    \begin{minipage}{0.48\textwidth}
        \centering
        \includegraphics[width=\textwidth]{spiky.png}
        \caption{Spiky kernel}
        \label{fig:spiky_mini}
    \end{minipage}
    \hfill
    \begin{minipage}{0.48\textwidth}
        \centering
        \includegraphics[width=\textwidth]{poly.png}
        \caption{Poly6 kernel}
        \label{fig:poly_mini}
    \end{minipage}
\end{figure}

\clearpage

## Surface Rendering via Ray Marching

### The Rendering Problem

The SPH produces a set of discrete particles, each representing a point within the fluid volume. A native rendering of these particles as individual spheres or points fails to convey the impression of continuous fluid with a coherent surface. Therefore, a method is required to reconstruct and render a visually plausible surface from this particle data in real-time.

### The Different techniques of Rendering Methods

Two primary techniques exist for this task: polygonization and screen-space rendering.

* Polygonization methods, such as the popular Marching Cubes algorithm, first convert the particle data into a density field on a 3D grid (a process known as voxelization). The algorithm then traverses this grid to generate a triangle mesh representing the fluid's surface, which can be rendered with standard techniques. While capable of producing high-quuality meshes, this approach was deemed unsuitable for this project due to tw main drawbacks: the high computational cost of voxelization and mesh generation, and the performance bottleneck associated with transferring a potentially large, dynamic mesh from the GPU(where the simulation runs) to the CPU for rendering in Unity.

* Screen-space methods, by contrast, operate directly on the GPU for each pixel on the screen. The choosen method for this project is Ray Marching an implicit surface, which avoids the creation of any intermediate geometry and is exceptionally well-suited for the parallel architecture of the GPU.

### Theory of Ray Marching

The core of the rendering technique is to define the fluid surface not with triangles, but as an implicit surface. The surface is defined as the set of all points $p$ in space where a function $\textbf{F(p)}$ equals zero. Specifically, we use a Signed Distance Field (SDF), a special type of implicit function that, for any point $p$, returns the shortest distance to the surface. The sign of the distance is positive outside the fluid and negative inside.

### Constructing the Fluid SDF

The SDF for the entire fluid is built by combining the SDFs of the individual particles. A single is presented by a sphere, whose SDF is given by:

equation here

To combine the fields from all particles into a single, smooth "blobby" surface, a simple minimum operation is insufficient as it would create sharp creases. Instead, a polynomial smooth minimum( $smin$ ) function is used. This smoothly bends the distance fields of nearby particles.

### The Ray Marching Algorithm

With the SDF defined, the surface can be rendered using a ray marching algorithm, also known as sphere tracing. For each pixel on the screen, a ray is cast from the camera. The algorithm proceeds iteratively:

1. From the ray's current position, the SDF is evaluated to find the distance $d$ to the surface.
2. The SDF guarantees that we can safely "march" the ray forward along its direction by the distance $d$ without passing through the surface.
3. This process is repeated. If $d$ becomes smaller than a small threshold (epsilon), the ray has hit the surface. If the ray travels too far, it is considered a miss.

This process is highly efficient as it takes the largest possible safe steps through empty space.

## Implementation and Methods

### Project Setup

A simple object scene was prepared in order for the simulation to be tested. There is a simple Ground cube and a few test objects which are the sphere and cube here.

We then create a SPH.cs script that initialize the particles that we are going to use.

### SPH.cs Script implementation

The First step in the script is to define a data structure that represents a single particle in the simulation. Each particle needs to store its physical properties (pressure, density), its current motion (velocity, force), and its position in the world. To make sure this structure is compatible between C# and the compute shader, it is laid out in memory sequentially. 

```csharp
[System.Serializable]
[StructLayout(LayoutKind.Sequential, Size=44)]
public struct Particle {
    public float pressure; // 4 bytes
    public float density; // 8 bytes
    public Vector3 currentForce; // 20 bytes
    public Vector3 velocity; // 32 bytes
    public Vector3 position; // 44 total bytes
}
```

This structure is only 44 bytes, which is small enough to handle thousands of particles effciently on the GPU. Each particle is initialized with default values in the script, and later updated every frame by the compute shaders (for density, pressure, force and integration).

We will then create a SPH class which will contain all the settings for the SPH Simulation

<ul>
<li>General Settings
    <ul>
        <li>collisionSphere:- We can assign our test sphere which can interact with the fluid.</li>
        <li>showSpheres:- Toggle to show the particles in the scene. (useful for debugging).</li>
        <li>numToSpawn:- Number of particles to spawn along each axis. (X,Y,Z).</li>
        <li>totalParticles:- The number of particles to spawn, which is set to each axis and multiply them.</li>
        <li>boxSize:- This is the bounding box of the scene.</li>
        <li>spawnCenter:- This is the spawn center of the particle grid in the scene.</li>
        <li>particleRadius:- This is radius of each of our particles.</li>
        <li>spawnJitter:- A value that will add randomness to our grid so that our grid of particles is not uniform.</li>
    </ul>
</li>
<li>Particle Rendering Settings
    <ul>
        <li>particleMesh:- This is the mesh of the particle which we will use.</li>
        <li>particleRenderSize:- The size in which the particle will render.</li>
        <li>material:- The material for the particles</li>
    </ul>
</li>
<li>Compute Settings
    <ul>
        <li>shader:- This is our compute shader.</li>
        <li>particles:- The particles list which stores all the particles in the scene.</li>
    </ul>
</li>
<li>Fluid Constants
    <ul>
        <li>boundDamping:- The value to use when the particles collide with the boundary, it reflects the velocity and scales it down to simulate energy loss</li>
        <li>viscosity:- The viscosity of our fluid.</li>
        <li>particleMass:- The individual mass of the particle.</li>
        <li>gasConstant:- The value of the gas constant.</li>
        <li>restDensity:- The rest density of the fluid.</li>
        <li>timestep:- The timestep of the scene.</li>
    </ul>
</li>
<li>Private Compute Shader Buffer Variables
    <ul>
        <li>_argsBuffer:- The arguments related to the values of the objects in the scene which we are sending to the compute shader.</li>
        <li>_particlesBuffer:- Sending the total amount of particles to the compute shader.</li>
        <li>integrateKernel:- This is the integrateKernel function from the compute shader.</li>
        <li>computeKernel:- This is the computeKernel function from the compute shader.</li>
        <li>densityKernel:- This is the densityKernel function from the compute shader.</li>
    </ul>
</li>
</ul>

```csharp
public class SPH : MonoBehaviour
{
    [Header("General")]
    public Transform collisionSphere;
    public bool showSpheres = true;
    public Vector3Int numToSpawn = new Vector3Int(10,10,10);
    private int totalParticles {
        get {
            return numToSpawn.x*numToSpawn.y*numToSpawn.z;
        }
    }
    public Vector3 boxSize = new Vector3(4,10,3);
    public Vector3 spawnCenter;
    public float particleRadius = 0.1f;
    public float spawnJitter = 0.2f;

    [Header("Particle Rendering")]
    public Mesh particleMesh;
    public float particleRenderSize = 8f;
    public Material material;
    
    [Header("Compute")]
    public ComputeShader shader;
    public Particle[] particles;

    [Header("Fluid Constants")]
    public float boundDamping = -0.3f;
    public float viscosity = -0.003f;
    public float particleMass = 1f;
    public float gasConstant = 2f;
    public float restingDensity = 1f;
    public float timestep = 0.007f;
    
    // Private Variables
    private ComputeBuffer _argsBuffer;
    public ComputeBuffer _particlesBuffer;
    private int integrateKernel;
    private int computeKernel;
    private int densityPressureKernel;
```