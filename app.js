'use strict';

var SwaggerExpress = require('swagger-express-mw');
var express = require('express');
const passport = require('passport');
var cors = require('cors');

var app = express();
app.use(cors());
app.use(require('express-session')({secret: process.env.SESSION_SECRET, resave: true, saveUninitialized: true }));
app.use(passport.initialize());
app.use(passport.session());

app.get('/auth/google', passport.authenticate('google', { scope: ['profile', 'email'] }));
app.get('/auth/google/callback', passport.authenticate('google', { failureRedirect: '/' }),
  (req, res) => {
    res.redirect('/');
  });
app.get('/logout', (req, res) => {
  req.logout();
  res.redirect('/');
});

function isAuthenticated(req, res, next) {
  if (req.isAuthenticated()) {
    return next();
  }
  res.status(401).send('Unauthorized');
}

app.get('/protected-api', isAuthenticated, (req, res) => {
  res.send('This is protected data');
});

module.exports = app;

var config = {
  appRoot: __dirname // required config
};

SwaggerExpress.create(config, function (err, swaggerExpress) {
  if (err) { throw err; }

  var basePath = swaggerExpress.runner.swagger.basePath;

  app.get('/', function(req, res) {
    res.redirect(basePath);
  });

  // define routes for documentation
  app.get(basePath, function (req, res, next) { // redirect to /docs with the correct schema
    res.redirect(basePath+'/docs?url='+basePath+'/swagger');
  });

  app.use(basePath+'/docs', express.static('node_modules/swagger-ui/dist'));

  // install swagger server middleware
  swaggerExpress.register(app);

  // redirect unimplemented routes to ebi atlas
  app.all(`${basePath}/*`, (req, res) => {
    // Construct the external URL, preserving the original path and query
    const ebiBase = 'https://www.ebi.ac.uk';
    const gxaUrl = ebiBase + req.originalUrl.replace(basePath,'');
    
    console.log(`Redirecting unhandled route ${req.originalUrl} to ${gxaUrl}`);
    res.redirect(301, gxaUrl); // 301 for permanent, 302 for temporary redirect
  });
  // start it up
  var version = +basePath.match(/\d+/);
  var port = 50003; // process.env.PORT || 10000 + version;
  app.listen(port);

  console.log('Listening on', port);
});

