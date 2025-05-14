cd /Users/akelbert/Developer/EMTF-FCU/database/xml;
outdir = '/Users/akelbert/Surveys/USGS/GIC/betaFactors';
tf = mttf.read;
tf.apresplt([1 1e5]); print(gcf,'-dpng',[outdir '/' tf.tfname '_apres.png']);
tf.impplt; print(gcf,'-dpng',[outdir '/' tf.tfname '_impedance.png']);