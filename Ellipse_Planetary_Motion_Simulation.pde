import peasy.*;

PeasyCam camera;
PImage sunTexture, backgroundTexture;
PImage mercuryTexture, venusTexture, earthTexture, marsTexture;
PImage jupiterTexture, saturnTexture, uranusTexture, neptuneTexture;

float scaleFactor = 1000;
float translationFactor = 1e8;
float rotationSpeedFactor = 1;
float visualSpeedFactor = 1;
float timeStep = 36000;

Planet[] planets;

class Planet {
  PImage texture;
  PVector pos;
  float radius;
  ArrayList<PVector> trail = new ArrayList<PVector>();
  float axialTilt;
  float rotationPeriod;
  float rotationAngle = 0;

  float semiMajorAxis;
  float eccentricity;
  float orbitalPeriod;
  float orbitalAngle = 0;

  int sweepInterval = 10;
  int lastSweepFrame = 0;
  PVector lastSweepPos;
  float lastArea = 0;
  PVector[] lastTriangle = null;

  ArrayList<Float> areaOverTime = new ArrayList<Float>();
  ArrayList<Float> velocityOverTime = new ArrayList<Float>();
  ArrayList<Float> distanceOverTime = new ArrayList<Float>();
  ArrayList<Float> timeStamps = new ArrayList<Float>();
  float totalTime = 0;

  Planet(PImage texture, float a, float e, float radiusKM, float tiltDeg, float siderealDaySec) {
    this.texture = texture;
    this.radius = radiusKM;
    this.axialTilt = tiltDeg;
    this.rotationPeriod = abs(siderealDaySec);
    this.semiMajorAxis = a;
    this.eccentricity = e;
    this.orbitalPeriod = TWO_PI * sqrt(pow(a, 3) / (6.67430e-11 * 1.989e30));
    this.pos = new PVector();
    this.lastSweepPos = new PVector();
  }

  void update() {
    totalTime += timeStep * visualSpeedFactor;

    float n = TWO_PI / orbitalPeriod;
    float M = (n * timeStep * visualSpeedFactor + orbitalAngle) % TWO_PI;

    float e = eccentricity;
    float E = M;
    for (int i = 0; i < 1; i++) {
      E = M + e * sin(E);
    }

    float theta = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
    
    float a = semiMajorAxis;
    float r = a * (1 - e * cos(E));
    pos = new PVector(r * cos(theta), r * sin(theta), 0);
    trail.add(pos.copy());

    orbitalAngle = M;
    rotationAngle += TWO_PI * timeStep / rotationPeriod * visualSpeedFactor;

    if (frameCount - lastSweepFrame >= sweepInterval) {
      float area = abs(lastSweepPos.x * pos.y - lastSweepPos.y * pos.x) / 2.0;
      lastArea = area;
      lastTriangle = new PVector[] {
        new PVector(0, 0, 0),
        lastSweepPos.copy(),
        pos.copy()
      };
      lastSweepPos = pos.copy();
      lastSweepFrame = frameCount;

      if (this == planets[2]) {
        float velocity = sqrt(6.67430e-11 * 1.989e30 * (2 / pos.mag() - 1 / semiMajorAxis));
        if (areaOverTime.size() > 0 && timeStamps.size() > 0) {
          float lastTime = timeStamps.get(timeStamps.size() - 1);
          if (totalTime - lastTime >= sweepInterval * timeStep * visualSpeedFactor) {
            areaOverTime.add(area);
            velocityOverTime.add(velocity);
            distanceOverTime.add(pos.mag());
            timeStamps.add(totalTime);
          }
        } else {
          areaOverTime.add(area);
          velocityOverTime.add(velocity);
          distanceOverTime.add(pos.mag());
          timeStamps.add(totalTime);
        }

        if (areaOverTime.size() > 300) {
          areaOverTime.remove(0);
          velocityOverTime.remove(0);
          distanceOverTime.remove(0);
          timeStamps.remove(0);
        }
      }
    }
  }

  void draw() {
    drawSphere(pos.x, pos.y, pos.z, texture, radius, axialTilt, rotationAngle);
    stroke(0, 150, 255);
    noFill();
    beginShape();
    for (PVector p : trail) {
      vertex(p.x / translationFactor, p.y / translationFactor, p.z / translationFactor);
    }
    endShape();
  }

  void drawSweptArea() {
    if (lastTriangle == null) return;
    PVector p0 = lastTriangle[0].copy().div(translationFactor);
    PVector p1 = lastTriangle[1].copy().div(translationFactor);
    PVector p2 = lastTriangle[2].copy().div(translationFactor);

    fill(255, 100, 0, 100);
    stroke(255, 150, 0);
    beginShape();
    vertex(p0.x, p0.y, p0.z);
    vertex(p1.x, p1.y, p1.z);
    vertex(p2.x, p2.y, p2.z);
    endShape(CLOSE);

    PVector centroid = new PVector((p0.x + p1.x + p2.x) / 3, (p0.y + p1.y + p2.y) / 3, (p0.z + p1.z + p2.z) / 3);
    pushMatrix();
    translate(centroid.x, centroid.y, centroid.z);
    fill(255);
    textSize(10);
    textAlign(CENTER);
    text(nf(lastArea / 1e16, 1, 2) + " ×10¹⁶ km²", 0, 0);
    popMatrix();
  }

  void drawInfo() {
    fill(255);
    textSize(10);
    PVector screenPos = pos.copy().div(translationFactor);
    pushMatrix();
    translate(screenPos.x, screenPos.y, screenPos.z);
    textAlign(LEFT);
    float velocity = sqrt(6.67430e-11 * 1.989e30 * (2 / pos.mag() - 1 / semiMajorAxis));
    text(
      "Orb. Period: " + nf(orbitalPeriod / (3600 * 24), 1, 1) + " days\n" +
      "Orb. Velocity: " + nf(velocity / 1000, 1, 2) + " km/s\n" +
      "Dist. from Sun: " + nf(pos.mag() / 1e9, 1, 0) + " million km",
      10, 0
    );
    popMatrix();
  }
}

void setup() {
  size(640, 360, P3D);
  sunTexture = loadImage("sunTexture.jpg");
  backgroundTexture = loadImage("stars.jpg");

  mercuryTexture = loadImage("mercuryTexture.jpg");
  venusTexture   = loadImage("venusTexture.jpg");
  earthTexture   = loadImage("earthTexture.jpg");
  marsTexture    = loadImage("marsTexture.jpg");
  jupiterTexture = loadImage("jupiterTexture.jpg");
  saturnTexture  = loadImage("saturnTexture.jpg");
  uranusTexture  = loadImage("uranusTexture.jpg");
  neptuneTexture = loadImage("neptuneTexture.jpg");

  float fov = PI / 3.0;
  float cameraZ = (height / 2.0) / tan(fov / 2.0);
  perspective(fov, float(width)/float(height), cameraZ / 10.0, cameraZ * 10000.0);
  camera = new PeasyCam(this, 0, 0, 0, 200);

  planets = new Planet[] {
    new Planet(mercuryTexture, 57.91e9, 0.2056, 2440, 0.03, 5070000),
    new Planet(venusTexture,   108.2e9, 0.0067, 3760, 177.4, -20997000),
    new Planet(earthTexture,   149.6e9, 0.0167, 6370, 23.44, 86164),
    new Planet(marsTexture,    227.9e9, 0.0934, 3390, 25.19, 88642),
    new Planet(jupiterTexture, 778.5e9, 0.0489, 71490, 3.13, 35730),
    new Planet(saturnTexture,  1.433e12, 0.0565, 58232, 26.73, 38362),
    new Planet(uranusTexture,  2.877e12, 0.0457, 25560, 97.77, -62064),
    new Planet(neptuneTexture, 4.503e12, 0.0113, 24764, 28.32, 57996)
  };
}

void draw() {
  background(0);
  drawSphere(0, 0, 0, sunTexture, 696340 / 2, 7.25, (frameCount * TWO_PI / 2192832) * 0.1);

  for (Planet p : planets) {
    p.update();
    p.draw();
    p.drawSweptArea();
    p.drawInfo();
  }

  drawGraphs(planets[2]); // Earth
}

void drawGraphs(Planet earth) {
  camera.beginHUD();

  int graphW = 200;
  int graphH = 50;
  int x0 = 10;
  int y0 = height - 3 * graphH - 60;

  // Compute dynamic scaling
  float maxArea = 0, maxVelocity = 0, maxDistance = 0, maxTime = 0;
  for (int i = 0; i < earth.areaOverTime.size(); i++) {
    maxArea = max(maxArea, earth.areaOverTime.get(i));
    maxVelocity = max(maxVelocity, earth.velocityOverTime.get(i));
    maxDistance = max(maxDistance, earth.distanceOverTime.get(i));
    maxTime = max(maxTime, earth.timeStamps.get(i));
  }
  
  maxArea *= 1.05;
  maxVelocity *= 1.05;
  maxDistance *= 1.05;
  maxTime *= 1.05;

  textAlign(LEFT);
  textSize(10);
  fill(255);

  // === Area vs Time ===
  fill(0);
  stroke(255);
  rect(x0, y0, graphW, graphH);
  text("Area vs Time", x0 + 5, y0 + 12);
  if (earth.areaOverTime.size() > 0) {
    text("Latest: " + nf(earth.areaOverTime.get(earth.areaOverTime.size() - 1) / 1e16, 1, 2) + " ×10¹⁶ km²", x0 + 5, y0 + 26);
  }
  
  
  noFill();
  stroke(255, 100, 0);
  beginShape();
  for (int i = 0; i < earth.areaOverTime.size(); i++) {
    float t = earth.timeStamps.get(i) / maxTime;
    float a = earth.areaOverTime.get(i) / maxArea;
    float px = x0 + t * graphW;
    float py = y0 + graphH - a * graphH;
    vertex(px, py);
  }
  endShape();

  // === Velocity vs Time ===
  y0 += graphH + 20;
  fill(0);
  stroke(255);
  rect(x0, y0, graphW, graphH);
  text("Velocity vs Time", x0 + 5, y0 + 12);
  if (earth.velocityOverTime.size() > 0) {
    text("Latest: " + nf(earth.velocityOverTime.get(earth.velocityOverTime.size() - 1) / 1000, 1, 2) + " km/s", x0 + 5, y0 + 26);
  }

  noFill();
  stroke(0, 200, 255);
  beginShape();
  for (int i = 0; i < earth.velocityOverTime.size(); i++) {
    float t = earth.timeStamps.get(i) / maxTime;
    float v = earth.velocityOverTime.get(i) / maxVelocity;
    float px = x0 + t * graphW;
    float py = y0 + graphH - v * graphH;
    vertex(px, py);
  }
  endShape();

  // === Distance vs Time ===
  y0 += graphH + 20;
  fill(0);
  stroke(255);
  rect(x0, y0, graphW, graphH);
  text("Distance vs Time", x0 + 5, y0 + 12);
  if (earth.distanceOverTime.size() > 0) {
    text("Latest: " + nf(earth.distanceOverTime.get(earth.distanceOverTime.size() - 1) / 1e9, 1, 2) + " million km", x0 + 5, y0 + 26);
  }

  noFill();
  stroke(0, 255, 100);
  beginShape();
  for (int i = 0; i < earth.distanceOverTime.size(); i++) {
    float t = earth.timeStamps.get(i) / maxTime;
    float d = earth.distanceOverTime.get(i) / maxDistance;
    float px = x0 + t * graphW;
    float py = y0 + graphH - d * graphH;
    vertex(px, py);
  }
  endShape();

  camera.endHUD();
}





void drawSphere(float x, float y, float z, PImage texture, float radiusKM, float axialTiltDeg, float rotationAngle) {
  x = x / translationFactor;
  y = y / translationFactor;
  z = z / translationFactor;

  noStroke();
  pushMatrix();
  translate(x, y, z);
  rotateZ(radians(axialTiltDeg));
  rotateY(rotationAngle * rotationSpeedFactor);

  textureMode(NORMAL);
  beginShape(TRIANGLE_STRIP);
  texture(texture);

  int detail = 60;
  float radius = radiusKM / scaleFactor;

  for (int i = 0; i <= detail; i++) {
    float theta1 = map(i, 0, detail, 0, PI);
    float theta2 = map(i + 1, 0, detail, 0, PI);

    for (int j = 0; j <= detail; j++) {
      float phi = map(j, 0, detail, 0, TWO_PI);

      float x1 = radius * sin(theta1) * cos(phi);
      float y1 = radius * cos(theta1);
      float z1 = radius * sin(theta1) * sin(phi);
      float u = map(j, 0, detail, 0, 1);
      float v1 = map(i, 0, detail, 0, 1);
      vertex(x1, y1, z1, u, v1);

      float x2 = radius * sin(theta2) * cos(phi);
      float y2 = radius * cos(theta2);
      float z2 = radius * sin(theta2) * sin(phi);
      float v2 = map(i + 1, 0, detail, 0, 1);
      vertex(x2, y2, z2, u, v2);
    }
  }

  endShape();
  popMatrix();
}
