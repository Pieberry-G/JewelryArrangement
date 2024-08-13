// HEAT METHOD

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

//#include "geometrycentral/surface/intrinsic_triangulation.h"
//#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
//#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "imgui.h"

#include "colormap.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr, mesh_uptr_diamond;
std::unique_ptr<VertexPositionGeometry> geometry_uptr, geometry_uptr_diamond;
// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh, *mesh_diamond;
VertexPositionGeometry* geometry, *geometry_diamond;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh, *psMesh_diamond;
std::string MESHNAME, MESHNAME_diamond;

// Some global variables
Vector<double> DELTA;                      // sources
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex
polyscope::SurfaceVertexColorQuantity* solnColors;
polyscope::SurfaceGraphQuantity* strokelines;
double maxPhi = 0.0;
double vertexRadius;
double isolinesRadius;
// Added by cyh
double strokelinesRadius;


void flipZ() {
    // Rotate mesh 180 deg about up-axis on startup
    glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
    for (Vertex v : mesh->vertices()) {
        Vector3 vec = geometry->inputVertexPositions[v];
        glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
        rvec = rot * rvec;
        geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
    }
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
}

/*
 * Map solution to mesh colors.
 */
std::vector<std::array<double, 3>> computeColors(const Vector<double>& sol) {

    // Determine maximum-magnitude element for scaling purposes
    maxPhi = 0;
    for (size_t i = 0; i < mesh->nVertices(); i++) {
        maxPhi = std::max(maxPhi, sol[i]); // distances should be >= 0
    }

    std::vector<std::array<double, 3>> colors;
    for (Vertex v : mesh->vertices()) {
        colors.push_back(mapToColor(maxPhi - sol[v.getIndex()], 0, maxPhi, "hot")); // invert colormap
    }
    return colors;
}

// Set mesh color.
void setColors(const std::vector<std::array<double, 3>>& colors) {
    solnColors = psMesh->addVertexColorQuantity("Solution", colors);
    solnColors->setEnabled(true);
}

// Added by cyh
/*
 * Display stroke.
 */
void showStroke() {

    std::vector<Vector3> positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    std::vector<glm::vec4>& stroke = polyscope::state::stroke;
    for (size_t i = 0; i < stroke.size(); i++) {
        positions.push_back({stroke[i].x, stroke[i].y, stroke[i].z});
        if (i > 0) {
            edgeInds.push_back({i - 1, i});
        }
    }
    strokelines = psMesh->addSurfaceGraphQuantity("Stroke", positions, edgeInds);
    strokelines->setEnabled(true);
    strokelines->setRadius(strokelinesRadius);
    strokelines->setColor({0.0, 0.0, 0.0});
}

/*
 * Show selected vertices.
 * This function gets called every time an element is selected on-screen.
 */
void showSelected() {

    // Show selected vertices in yellow
    std::vector<Vector3> vertPos;
    std::vector<std::array<size_t, 2>> vertInd;
    DELTA = Vector<double>::Zero(mesh->nVertices());
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
         it != polyscope::state::subset.vertices.end(); ++it) {
        vertPos.push_back(geometry->inputVertexPositions[*it]);
        DELTA[*it] = 1;
    }
    polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    showVerts->setEnabled(true);
    showVerts->setRadius(vertexRadius);
    showVerts->setColor({1.0, 0.65, 0.0});

    // Show the currently selected vertex in red.
    int currIdx = polyscope::state::currVertexIndex;
    if (currIdx != -1) {
        std::vector<Vector3> pos = {geometry->inputVertexPositions[currIdx]};
        currVert = psMesh->addSurfaceGraphQuantity("current vertex", pos, std::vector<std::array<size_t, 2>>());
        currVert->setEnabled(true);
        currVert->setRadius(vertexRadius);
        currVert->setColor({1.0, 0.0, 0.0});
    } else {
        currVert->setEnabled(false);
    }
}

void redraw() {
    showSelected();
    polyscope::requestRedraw();
}

void functionCallback() {

    if (ImGui::Button("Reset")) {
        polyscope::state::subset.vertices.clear();
        polyscope::state::currVertexIndex = -1;
        polyscope::state::stroke.clear();
        showStroke();
        psMesh->setSurfaceColor({1.0, 0.45, 0.0});
        solnColors->setEnabled(false);
        redraw();
    }
}


int main() {

    // If a mesh name was not given, use default mesh.
    std::string filepath;
    //filepath = "../input/bunny.obj";
    filepath = "../input/rings/ring1.obj";
    filepath = "../input/rings/ring_18.obj";

    MESHNAME = polyscope::guessNiceNameFromPath(filepath);

    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    // Get indices for element picking
    polyscope::state::facePickIndStart = mesh->nVertices();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    // Initalize
    flipZ();
    double lengthScale = geometry->meanEdgeLength();
    polyscope::state::edgeLengthScale = lengthScale;
    vertexRadius = lengthScale * 0.2;
    strokelinesRadius = 0.005f;
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange
    currVert =
        psMesh->addSurfaceGraphQuantity("current vertex", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    // initialize to something in case "Reset" is pressed before anything happens
    solnColors = psMesh->addVertexColorQuantity("Solution", std::vector<std::array<double, 3>>(mesh->nVertices()));
    DELTA = Vector<double>::Zero(mesh->nVertices());








    // Load diamond mesh.
    std::string diamondFilepath = "../input/diamonds/diamond3.obj";

    MESHNAME_diamond = polyscope::guessNiceNameFromPath(diamondFilepath);

    // Load mesh
    std::tie(mesh_uptr_diamond, geometry_uptr_diamond) = readManifoldSurfaceMesh(diamondFilepath);
    mesh_diamond = mesh_uptr_diamond.release();
    geometry_diamond = geometry_uptr_diamond.release();

    psMesh_diamond = polyscope::registerSurfaceMesh(MESHNAME_diamond, geometry_diamond->inputVertexPositions, mesh_diamond->getFaceVertexList(),
                                       polyscopePermutations(*mesh_diamond));

    // Give control to the polyscope gui
    polyscope::show();

    delete mesh;
    delete geometry;

    delete mesh_diamond;
    delete geometry_diamond;

    return EXIT_SUCCESS;
}