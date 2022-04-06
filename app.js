'use strict';

var SwaggerExpress = require('swagger-express-mw');
var express = require('express');
var cors = require('cors');

var app = express();
app.use(cors());

//var fs = require('fs');
//var yaml = require('js-yaml');

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

  // start it up
  var version = +basePath.match(/\d+/);
  var port = 20003; // process.env.PORT || 10000 + version;
  app.listen(port);

  console.log('Listening on', port);
});

