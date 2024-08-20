#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

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
double edgeRadius;
// Added by cyh
double strokelinesRadius;

glm::vec<3, float> BLUE_VEC = { 0.11, 0.388, 0.89 };
std::array<double, 3> BLUE = { 0.11, 0.388, 0.89 };
glm::vec<3, float> ORANGE_VEC = { 1.0, 0.45, 0.0 };
std::array<double, 3> ORANGE = { 1.0, 0.45, 0.0 };


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

    //// Show selected vertices.
    //std::vector<Vector3> vertPos;
    //std::vector<std::array<size_t, 2>> vertInd;
    //for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
    //     it != polyscope::state::subset.vertices.end(); ++it) {
    //    vertPos.push_back(geometry->inputVertexPositions[*it]);
    //}
    //polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    //showVerts->setEnabled(true);
    //showVerts->setRadius(vertexRadius);
    //showVerts->setColor(ORANGE_VEC);

    //// Show selected edges.
    //std::vector<Vector3> edgePos;
    //std::vector<std::array<size_t, 2>> edgeInd;
    //for (std::set<size_t>::iterator it = polyscope::state::subset.edges.begin();
    //     it != polyscope::state::subset.edges.end(); ++it) {
    //    Edge e = mesh->edge(*it);
    //    edgePos.push_back(geometry->inputVertexPositions[e.firstVertex()]);
    //    edgePos.push_back(geometry->inputVertexPositions[e.secondVertex()]);
    //    size_t i = edgeInd.size();
    //    edgeInd.push_back({2 * i, 2 * i + 1});
    //}
    //polyscope::SurfaceGraphQuantity* showEdges = psMesh->addSurfaceGraphQuantity("selected edges", edgePos, edgeInd);
    //showEdges->setEnabled(true);
    //showEdges->setRadius(edgeRadius);
    //showEdges->setColor(BLUE_VEC);

    // Show selected faces.
    std::vector<std::array<double, 3>> faceColors(mesh->nFaces());
    for (size_t i = 0; i < mesh->nFaces(); i++) {
        faceColors[i] = ORANGE;
    }
    for (std::set<size_t>::iterator it = polyscope::state::subset.faces.begin();
         it != polyscope::state::subset.faces.end(); ++it) {
        faceColors[*it] = {0.5, 0, 0.5};
    }
    polyscope::SurfaceFaceColorQuantity* showFaces = psMesh->addFaceColorQuantity("selected faces", faceColors);
    showFaces->setEnabled(true);
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
    vertexRadius = 0.005f;
    vertexRadius = 0.005f;
    strokelinesRadius = 0.005f;
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor(ORANGE_VEC); // orange
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