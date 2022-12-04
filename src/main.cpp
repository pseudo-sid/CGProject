#include "utils.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <queue>

#define GLM_FORCE_RADIANS
#define GLM_ENABLE_EXPERIMENTAL

#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

std::vector<std::vector<int>> crownPoints {{-101, 27}, {-102, 27}, {-103, 26}, {-104, 26}, {-105, 26}, {-106, 26}, {-107, 25}, {-108, 25}, {-109, 24}, {-109, 23}, {-110, 22}, {-111, 22}, {-111, 21}, {-111, 20}, {-111, 19}, {-111, 18}, {-112, 18}, {-112, 17}, {-112, 16}, {-112, 15}, {-112, 14}, {-112, 13}, {-112, 12}, {-112, 11}, {-111, 10}, {-111, 9}, {-111, 7}, {-110, 5}, {-110, 4}, {-110, 2}, {-110, 0}, {-110, -1}, {-110, -3}, {-110, -5}, {-110, -6}, {-110, -7}, {-110, -8}, {-110, -10}, {-110, -11}, {-111, -11}, {-112, -12}, {-112, -14}, {-113, -14}, {-113, -15}, {-114, -16}, {-114, -17}, {-114, -18}, {-115, -18}, {-117, -19}, {-117, -20}, {-118, -21}, {-118, -22}, {-120, -22}, {-120, -23}, {-121, -24}, {-122, -25}, {-122, -26}, {-123, -27}, {-124, -28}, {-124, -29}, {-124, -30}, {-124, -31}, {-124, -32}, {-124, -33}, {-123, -34}, {-122, -36}, {-122, -38}, {-121, -39}, {-120, -40}, {-119, -42}, {-119, -44}, {-118, -45}, {-118, -46}, {-117, -47}, {-115, -49}, {-115, -50}, {-114, -50}, {-114, -52}, {-113, -54}, {-111, -55}, {-111, -56}, {-110, -58}, {-109, -58}, {-108, -60}, {-107, -60}, {-106, -61}, {-106, -62}, {-105, -62}, {-104, -62}, {-103, -63}, {-103, -64}, {-102, -65}, {-102, -66}, {-102, -67}, {-102, -68}, {-101, -68}, {-101, -69}, {-100, -69}, {-100, -70}, {-99, -70}, {-99, -71}, {-98, -72}, {-98, -73}, {-97, -74}, {-95, -74}, {-94, -75}, {-93, -75}, {-92, -76}, {-90, -77}, {-89, -77}, {-88, -77}, {-86, -78}, {-84, -78}, {-83, -78}, {-82, -78}, {-80, -78}, {-79, -78}, {-78, -78}, {-76, -78}, {-75, -78}, {-74, -78}, {-73, -78}, {-71, -78}, {-70, -78}, {-70, -77}, {-69, -77}, {-67, -77}, {-66, -77}, {-64, -77}, {-63, -77}, {-62, -77}, {-60, -77}, {-59, -77}, {-58, -77}, {-56, -77}, {-55, -77}, {-54, -78}, {-53, -78}, {-52, -78}, {-52, -79}, {-51, -80}, {-50, -80}, {-50, -81}, {-49, -82}, {-48, -83}, {-47, -83}, {-47, -84}, {-46, -85}, {-46, -86}, {-45, -86}, {-44, -86}, {-43, -87}, {-42, -88}, {-41, -89}, {-40, -90}, {-39, -90}, {-38, -90}, {-36, -91}, {-34, -92}, {-33, -92}, {-31, -92}, {-30, -93}, {-29, -93}, {-27, -94}, {-26, -94}, {-24, -94}, {-22, -94}, {-21, -94}, {-19, -94}, {-18, -94}, {-17, -94}, {-15, -94}, {-14, -94}, {-13, -94}, {-11, -94}, {-10, -94}, {-9, -93}, {-8, -93}, {-6, -93}, {-5, -92}, {-2, -92}, {-1, -92}, {2, -92}, {5, -92}, {8, -92}, {10, -92}, {14, -92}, {17, -92}, {19, -92}, {22, -92}, {24, -92}, {26, -92}, {28, -92}, {30, -92}, {32, -93}, {33, -93}, {34, -93}, {37, -94}, {39, -94}, {41, -94}, {42, -94}, {44, -94}, {46, -95}, {47, -95}, {49, -95}, {50, -95}, {51, -95}, {52, -95}, {53, -94}, {54, -94}, {55, -93}, {55, -92}, {55, -91}, {56, -90}, {57, -90}, {58, -89}, {58, -88}, {26, 24}, {25, 24}, {25, 23}, {25, 22}, {24, 22}, {22, 22}, {22, 21}, {22, 20}, {21, 20}, {21, 19}, {20, 18}, {19, 18}, {18, 17}, {18, 16}, {17, 15}, {16, 15}, {16, 14}, {14, 14}, {13, 14}, {11, 13}, {10, 13}, {9, 12}, {6, 11}, {4, 11}, {2, 10}, {1, 10}, {0, 10}, {-1, 10}, {-2, 10}, {-3, 10}, {-5, 10}, {-6, 10}, {-7, 10}, {-8, 10}, {-9, 10}, {-10, 10}, {-10, 9}, {-11, 9}, {-11, 8}, {-12, 8}, {-13, 8}, {-14, 8}, {-14, 7}, {-15, 7}, {-16, 7}, {-17, 7}, {-18, 8}, {-19, 9}, {-19, 10}, {-20, 10}, {-22, 11}, {-22, 13}, {-23, 13}, {-25, 14}, {-26, 16}, {-27, 18}, {-28, 19}, {-30, 19}, {-30, 21}, {-31, 22}, {-32, 22}, {-33, 22}, {-34, 23}, {-34, 24}, {-35, 24}, {-36, 24}, {-37, 24}, {-38, 25}, {-39, 26}, {-40, 26}, {-41, 26}, {-42, 26}, {-42, 27}, {-42, 28}, {-44, 28}, {-45, 29}, {-46, 29}, {-47, 29}, {-48, 29}, {-49, 29}, {-49, 28}, {-50, 28}, {-50, 27}, {-51, 27}, {-52, 27}, {-52, 26}, {-53, 26}, {-54, 26}, {-55, 26}, {-57, 26}, {-58, 27}, {-59, 27}, {-60, 28}, {-62, 29}, {-64, 30}, {-66, 30}, {-67, 31}, {-69, 32}, {-71, 33}, {-72, 33}, {-74, 34}, {-76, 34}, {-77, 34}, {-78, 34}, {-79, 34}, {-80, 34}, {-81, 34}, {-82, 34}, {-83, 34}, {-84, 34}, {-85, 34}, {-86, 34}, {-87, 34}, {-88, 34}, {-88, 33}, {-89, 33}, {-90, 32}, {-90, 31}, {-91, 30}, {-92, 30}, {-93, 30}, {-94, 29}, {-94, 28}, {-94, 27}, {-95, 26}, {-96, 26}, {-97, 26}, {-98, 26}, {-99, 26}, {-100, 26}, {-101, 26}, {-102, 27}, {-102, 28}, {-103, 28}, {-104, 28}, {-104, 29}, {-105, 29}};
std::vector<std::vector<int>> trunkPoints {{-12, -10}, {-9, -2}, {-9, 6}, {-8, 21}, {-9, 34}, {-10, 42}, {-12, 49}, {-13, 58}, {-14, 65}, {-14, 70}, {-14, 82}, {-13, 88}, {-9, 97}, {-9, 101}, {-9, 109}};

int screen_width = 1000, screen_height = 1000;
GLint vModel_uniform, vView_uniform, vProjection_uniform; // GL_INT
GLint vColor_uniform; // GL_INT
glm::mat4 modelT, viewT, projectionT; // model, view and projection transformations

double oldX, oldY, currentX, currentY; // last mouse x position, last mouse y position, current mouse x position, current mouse y position
bool isDragging = false;

void setupModelTransformation(unsigned int &);
void setupViewTransformation(unsigned int &);
void setupProjectionTransformation(unsigned int &);
void createPoints(unsigned int &, unsigned int &, std::vector<std::vector<int>>&);
void createLine(unsigned int &, unsigned int &, const std::vector<int>&, const std::vector<int>&);
glm::vec3 getTrackBallVector(double x, double y);

std::vector<int> findMinMaxPositions(std::vector<std::vector<int>>& points) {
    int min_x = INT_MAX;
    int max_x = INT_MIN;
    int min_y = INT_MAX;
    int max_y = INT_MIN;

    for (auto p : points) {
        min_x = std::min(p[0], min_x);
        max_x = std::max(p[0], max_x);
        min_y = std::min(p[1], min_y);
        max_y = std::max(p[1], max_y);
    }

    std::vector<int> result { min_x, max_x, min_y, max_y };
    return result;
}

std::vector<std::vector<int>> pointsInsideCurve(std::vector<std::vector<int>>& curvePoints, std::vector<std::vector<int>> randomPoints) {
    std::vector<std::vector<int>> insidePoints;
    std::vector<int> maxMins = findMinMaxPositions(curvePoints);
    int min_x = maxMins[1];
    int max_x = maxMins[0];
    std::unordered_map<int, std::pair<int, int>> rangeY;
    
    // find minimum and maximum value of x, and consecutively find min and max values of Y for every x
    for (std::vector<int> &point: curvePoints) {
        int x = point[0], y = point[1];
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        
        // if x is not a key of hashmap, initialize values corresponding to x
        if (rangeY.find(x) == rangeY.end()) {
            rangeY[x].first = maxMins[1];
            rangeY[x].second = maxMins[0];
        }

        rangeY[x].first = std::min(rangeY[x].first, y);
        rangeY[x].second = std::max(rangeY[x].second, y);
    }

    for (std::vector<int> &point: randomPoints) {
        int x = point[0], y = point[1];

        if (x >= min_x and x <= max_x and y >= rangeY[x].first and y <= rangeY[x].second)
            insidePoints.push_back({x, y});
    }

    return insidePoints;
}

std::vector<std::vector<int>> generateRandomPoints(int min_x, int max_x, int min_y, int max_y, int n) {
    std::vector<std::vector<int>> points;

    for (int i = 0; i < n; i++) {
        int x = min_x + (std::rand() % (max_x - min_x + 1));
        int y = min_y + (std::rand() % (max_y - min_y + 1));
        points.push_back(std::vector<int>{x, y});
    }

    return points;
}

std::string randomString(const int length) {
    static const char characterPool[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::string tempString;
    tempString.reserve(length);

    for (int i = 0; i < length; ++i)
        tempString += characterPool[rand() % (sizeof(characterPool) - 1)];
    
    return tempString;
}

double distance(const std::vector<int>& a, const std::vector<int>& b) {
    return pow(
        pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2), 
        0.5
    );
}

std::unordered_map<std::string, double> distanceFromPoint(const std::vector<int>& start_point, std::unordered_set<std::string>& keys, std::unordered_map<std::string, std::vector<int>>& points) {
    std::unordered_map<std::string, double> distances;

    for (auto key : keys)
        distances[key] = distance(start_point, points[key]);

    return distances;
}

class Element {
public:
    std::vector<int> point;
    std::string key;

    Element(std::vector<int> p, std::string k) {
        point = p;
        key = k;
    }
};

void createLine(unsigned int & program, unsigned int & shape_VAO, const std::vector<int>& a, const std::vector<int>& b) {
    glUseProgram(program);

    // Bind shader variables
    int vVertex_attrib = glGetAttribLocation(program, "vVertex");
    
    if (vVertex_attrib == -1) {
        fprintf(stderr, "Could not bind location: vVertex\n");
        exit(0);
    }

    //Generate VAO object
    glGenVertexArrays(1, &shape_VAO); // creating a new Vertex Arrays Object
    glBindVertexArray(shape_VAO); // binding to the VAO

    // Create VBOs for the VAO
    // Position information (data + format)
    int nVertices = 2 * 3;
    GLfloat *expanded_vertices = new GLfloat[nVertices]; // GL_FLOAT
    int multiplier = 30;

    expanded_vertices[0] = (GLfloat) a[0] * 0.1;
    expanded_vertices[1] = (GLfloat) a[1] * 0.1;
    expanded_vertices[2] = (GLfloat) a[2] * 0.1;
    expanded_vertices[3] = (GLfloat) b[0] * 0.1;
    expanded_vertices[4] = (GLfloat) b[1] * 0.1;
    expanded_vertices[5] = (GLfloat) b[2] * 0.1;

    GLuint vertex_VBO;
    glGenBuffers(1, &vertex_VBO); // creating a vector buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO); // binding the vector buffer or selecting it
    glBufferData(GL_ARRAY_BUFFER, nVertices * sizeof(GLfloat), expanded_vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vVertex_attrib);
    glVertexAttribPointer(
        vVertex_attrib,
        3, // X, Y, Z
        GL_FLOAT,
        GL_FALSE,
        3 * sizeof(GLfloat),
        0
    );
    delete []expanded_vertices;

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0); //Unbind the VAO to disable changes outside this function.
}

std::unordered_map<std::string, std::string> buildTreeSkeleton(std::vector<int> _start, std::string _start_key, std::unordered_map<std::string, std::vector<int>> points) {
    std::unordered_set<std::string> remaining;
    std::unordered_set<std::string> taken;
    std::queue<Element*> q;
    std::unordered_map<std::string, std::string> parent;

    for (auto p : points)
        remaining.insert(p.first);

    q.push(new Element(_start, _start_key));

    for (auto s : remaining)
        parent[s] = s;

    parent[_start_key] = _start_key;
    
    while (!q.empty()) {
        Element* top = q.front();
        q.pop();

        std::unordered_map<std::string, double> distances = distanceFromPoint(top->point, remaining, points);
        std::vector<std::pair<std::string, double>> items;

        for (auto e : distances)
            items.push_back(e);

        sort(items.begin(), items.end(), [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
            return a.second < b.second;
        });
        int n = std::min((int) (3 + rand() % (5 - 3 + 1)), (int) items.size());
        items.resize(n);
        std::random_shuffle(items.begin(), items.end());
        
        for (int i = 0; i < n; i++) {
            taken.insert(items[i].first);
            remaining.erase(items[i].first);
        }

        for (int i = 0; i < n; i++) {
            if (taken.find(items[i].first) != taken.end()) {
                parent[items[i].first] = top->key;
                q.push(new Element(points[items[i].first], items[i].first));
            }
        }
    }

    // for (auto e : parent) {
    //     std::cout << e.first << " parent is " << e.second << std::endl;
    // }

    return parent;
}

std::pair<std::string, std::vector<int>> minPoint(std::unordered_map<std::string, std::vector<int>> points) {
    std::vector<std::pair<std::string, std::vector<int>>> items;
    std::pair<std::string, std::vector<int>> result;

    for (auto e : points)
        items.push_back(e);

    result.first = items[0].first;
    result.second = items[0].second;

    for (int i = 1; i < items.size(); i++) {
        if (items[i].second[1] < result.second[1]) {
            result.first = items[i].first;
            result.second = items[i].second;
        }
    }

    return result;
}

void addZCoordinateToPoints(std::vector<std::vector<int>>& points, int range) {
    for (int i = 0; i < points.size(); i++) {
        int z = (int) (0 + rand() % ((2 * range) + 1)) - range;
        points[i].push_back(z);
        // std::cout << z << std::endl;
    }
}

int main(int, char **) {
    std::vector<int> maxMins = findMinMaxPositions(crownPoints);
    int xLength = maxMins[1] - maxMins[0];
    int yLength = maxMins[3] - maxMins[2];
    std::vector<std::vector<int>> randomPoints = generateRandomPoints(maxMins[0], maxMins[1], maxMins[2], maxMins[3], 300);
    std::vector<std::vector<int>> insideCrownPoints = pointsInsideCurve(crownPoints, randomPoints);

    addZCoordinateToPoints(insideCrownPoints, xLength * 0.5);
    addZCoordinateToPoints(trunkPoints, 0);
    // addZCoordinateToPoints();

    std::unordered_map<std::string, std::vector<int>> crownPointsMap;
    std::unordered_map<std::string, std::vector<int>> trunkPointsMap;
    std::unordered_map<std::string, std::vector<int>> insideCrownPointsMap;

    for (auto p : crownPoints)
        crownPointsMap[randomString(25)] = p;

    for (auto p : trunkPoints)
        trunkPointsMap[randomString(25)] = p;

    for (auto p : insideCrownPoints)
        insideCrownPointsMap[randomString(25)] = p;

    std::pair<std::string, std::vector<int>> mp = minPoint(trunkPointsMap);
    std::unordered_map<std::string, std::string> branchParents = buildTreeSkeleton(mp.second, mp.first, insideCrownPointsMap);

    // setup window
    GLFWwindow *window = setupWindow(screen_width, screen_height);
    ImGuiIO &io = ImGui::GetIO(); // create IO object

    ImVec4 clearColor = ImVec4(1.0f, 1.0f, 1.0f, 1.00f); // black color

    unsigned int shaderProgram = createProgram("./shaders/vshader.vs", "./shaders/fshader.fs"); // loading vector and fragment shaders source codes, and create a program

    // get handle to color variable in shader
    vColor_uniform = glGetUniformLocation(shaderProgram, "vColor"); // location of a uniform variable vColor

    if (vColor_uniform == -1) {
        fprintf(stderr, "Could not bind location: vColor\n");
        exit(0);
    }

    glUseProgram(shaderProgram); // activating the shader program

    unsigned int VAO, VBO, VCO;
    unsigned int insideCrownVAOs [branchParents.size()];
    glGenVertexArrays(1, &VAO);
    createPoints(shaderProgram, VAO, crownPoints);
    glGenVertexArrays(1, &VBO);
    createPoints(shaderProgram, VBO, trunkPoints);

    int iii = 0;

    for (auto e : branchParents) {
        std::string child = e.first;
        std::string parent = e.second;
        std::vector<int> childPoint = child == mp.first ? mp.second : insideCrownPointsMap[child];
        std::vector<int> parentPoint = parent == mp.first ? mp.second : insideCrownPointsMap[parent];
        
        glGenVertexArrays(1, &insideCrownVAOs[iii]);
        createLine(shaderProgram, insideCrownVAOs[iii], childPoint, parentPoint);
        iii++;
    }

    setupModelTransformation(shaderProgram);
    setupViewTransformation(shaderProgram);
    setupProjectionTransformation(shaderProgram);

    oldX = oldY = currentX = currentY = 0.0;
    int prevLeftButtonState = GLFW_RELEASE;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Get current mouse position
        int leftButtonState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        double x, y;
        glfwGetCursorPos(window, &x, &y);

        if (leftButtonState == GLFW_PRESS && prevLeftButtonState == GLFW_RELEASE) {
            isDragging = true;
            currentX = oldX = x;
            currentY = oldY = y;
        } else if (leftButtonState == GLFW_PRESS && prevLeftButtonState == GLFW_PRESS) {
            currentX = x;
            currentY = y;
        } else if (leftButtonState == GLFW_RELEASE && prevLeftButtonState == GLFW_PRESS) {
            isDragging = false;
        }

        // Rotate based on mouse drag movement
        prevLeftButtonState = leftButtonState;
        
        if (isDragging && (currentX != oldX || currentY != oldY)) {
            glm::vec3 va = getTrackBallVector(oldX, oldY);
            glm::vec3 vb = getTrackBallVector(currentX, currentY);

            float angle = acos(std::min(1.0f, glm::dot(va, vb)));
            glm::vec3 axis_in_camera_coord = glm::cross(va, vb);
            glm::mat3 camera2object = glm::inverse(glm::mat3(viewT * modelT));
            glm::vec3 axis_in_object_coord = camera2object * axis_in_camera_coord;
            modelT = glm::rotate(modelT, angle, axis_in_object_coord);
            glUniformMatrix4fv(vModel_uniform, 1, GL_FALSE, glm::value_ptr(modelT));

            oldX = currentX;
            oldY = currentY;
        }

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        glUseProgram(shaderProgram);

        {
            ImGui::Begin("Information");
            ImGui::Text("%.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clearColor.x, clearColor.y, clearColor.z, clearColor.w);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glBindVertexArray(VAO);
        glUniform3f(vColor_uniform, 0.5, 0.5, 0.5);
        glPointSize(2);
        glDrawArrays(GL_POINTS, 0, crownPoints.size());

        glBindVertexArray(VBO);
        glUniform3f(vColor_uniform, 1.0, 0.5, 0.5);
        glPointSize(2);
        glDrawArrays(GL_POINTS, 0, trunkPoints.size());

        for (int i = 0; i < branchParents.size(); i++) {
            glBindVertexArray(insideCrownVAOs[i]);
            glUniform3f(vColor_uniform, 0., 0., 0.);
            glDrawArrays(GL_LINES, 0, 2);
        }

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    cleanup(window);
    return 0;
}

void createPoints(unsigned int &program, unsigned int &shape_VAO, std::vector<std::vector<int>>& points) {
    glUseProgram(program);

    // Bind shader variables
    int vVertex_attrib = glGetAttribLocation(program, "vVertex");
    
    if (vVertex_attrib == -1) {
        fprintf(stderr, "Could not bind location: vVertex\n");
        exit(0);
    }

    //Generate VAO object
    glGenVertexArrays(1, &shape_VAO); // creating a new Vertex Arrays Object
    glBindVertexArray(shape_VAO); // binding to the VAO

    // Create VBOs for the VAO
    // Position information (data + format)
    int nVertices = points.size() * 3;
    GLfloat *expanded_vertices = new GLfloat[nVertices]; // GL_FLOAT
    int multiplier = 30;

    for (int i = 0; i < points.size(); i++) {
        expanded_vertices[i * 3] = (GLfloat) points[i][0] * 0.1;
        expanded_vertices[i * 3 + 1] = (GLfloat) points[i][1] * 0.1;
        expanded_vertices[i * 3 + 2] = (GLfloat) points[i][2] * 0.1;
    }

    GLuint vertex_VBO;
    glGenBuffers(1, &vertex_VBO); // creating a vector buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO); // binding the vector buffer or selecting it
    glBufferData(GL_ARRAY_BUFFER, nVertices * sizeof(GLfloat), expanded_vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vVertex_attrib);
    glVertexAttribPointer(
        vVertex_attrib,
        3, // X, Y, Z
        GL_FLOAT,
        GL_FALSE,
        3 * sizeof(GLfloat),
        0
    );
    delete []expanded_vertices;

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0); //Unbind the VAO to disable changes outside this function.
}

void setupModelTransformation(unsigned int &program)
{
    // Modelling transformations (Model -> World coordinates)
    modelT = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, 0.0)); // Model coordinates are the world coordinates

    // Pass on the modelling matrix to the vertex shader
    glUseProgram(program);
    vModel_uniform = glGetUniformLocation(program, "vModel");

    if (vModel_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: vModel\n");
        exit(0);
    }

    glUniformMatrix4fv(vModel_uniform, 1, GL_FALSE, glm::value_ptr(modelT));
}

void setupViewTransformation(unsigned int &program)
{
    // Viewing transformations (World -> Camera coordinates
    // Camera at (0, 0, 100) looking down the negative Z-axis in a right handed coordinate system
    viewT = glm::lookAt(glm::vec3(0.0, 0.0, 40.0), glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 1.0, 0.0));

    // Pass-on the viewing matrix to the vertex shader
    glUseProgram(program);
    vView_uniform = glGetUniformLocation(program, "vView");
    if (vView_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: vView\n");
        exit(0);
    }
    glUniformMatrix4fv(vView_uniform, 1, GL_FALSE, glm::value_ptr(viewT));
}

void setupProjectionTransformation(unsigned int &program)
{
    // Projection transformation
    projectionT = glm::perspective(45.0f, (GLfloat)screen_width / (GLfloat)screen_height, 0.1f, 1000.0f);

    // Pass on the projection matrix to the vertex shader
    glUseProgram(program);
    vProjection_uniform = glGetUniformLocation(program, "vProjection");
    if (vProjection_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: vProjection\n");
        exit(0);
    }
    glUniformMatrix4fv(vProjection_uniform, 1, GL_FALSE, glm::value_ptr(projectionT));
}

glm::vec3 getTrackBallVector(double x, double y)
{
    glm::vec3 p = glm::vec3(2.0 * x / screen_width - 1.0, 2.0 * y / screen_height - 1.0, 0.0); // Normalize to [-1, +1]
    p.y = -p.y;                                                                                // Invert Y since screen coordinate and OpenGL coordinates have different Y directions.

    float mag2 = p.x * p.x + p.y * p.y;
    if (mag2 <= 1.0f)
        p.z = sqrtf(1.0f - mag2);
    else
        p = glm::normalize(p); // Nearest point, close to the sides of the trackball
    return p;
}
