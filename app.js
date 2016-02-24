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

  // install middleware
  swaggerExpress.register(app);

  var port = process.env.PORT || 10011;
  app.listen(port);

  console.log('Listening on', port);

  // define routes
  app.get('/', function (req, res, next) { // redirect to /docs with the correct schema
    res.redirect('/docs?url=/swagger');
  });

  app.use('/docs', express.static('node_modules/swagger-ui/dist'));
});

