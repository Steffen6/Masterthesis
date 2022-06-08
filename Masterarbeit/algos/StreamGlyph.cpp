#include <fantom/algorithm.hpp>
#include <fantom/graphics.hpp>
#include <fantom/register.hpp>
#include <fantom/dataset.hpp>
#include <fantom/math.hpp>
#include <fantom/cells.hpp>
#include<fantom/datastructures/ValueArray.hpp>
#include<fantom/datastructures/domains/LineSet.hpp>
#include <fantom/datastructures/interfaces/Field.hpp>
// needed for BoundinSphere-Calculation and normal calculation
//#include <fantom/graphics/HelperFunctions.hpp>
#include <utils/Graphics/include/HelperFunctions.hpp>
#include <utils/Graphics/include/ObjectRenderer.hpp>
#include <utils/math/include/eigenvalues.hpp>
#include <vector>
#include <stdexcept>
#include <math.h>

//source/toolboxes/utils/math/include/eigenvalues.hpp
//Mit GetDomain oder das Lineset bekommen.
//Als nächstes Rechteckige Röhren um die Linien bauen
//Dataalgorithmus machen und dann den ShowSurface-Algorithmus von Grid nutzen
//DomainFactory um das Dreiecksgitter zu erzeugen
//Farbgebung anschauen ColorMapping und Transferfunktion
//28.04 15:00 Uhr
using namespace fantom;
namespace{

    class Masterarbeit : public DataAlgorithm{
        public:

        struct Options : public DataAlgorithm::Options
        {
            Options( fantom::Options::Control& control )
                : DataAlgorithm::Options( control )
            {
                //Convert cellwise to pointwise mit direct, NICHT smooth
                add<Field<3, Matrix3> >("Field", "A 3D Matrix field", definedOn<Grid<3> >(Grid<3>::Points));//Tensor<double,3,3>
                add< LineSet< 3 > >( "LineSet","The Lineset of the StreamLine");
                add< double >( "Scale", "The Scale for the things", 0.1 );
                add< double >( "MinValue", "The Minimum Value for Eigenvectors", 1.0 );
            }
        };

        struct DataOutputs : public DataAlgorithm::DataOutputs {
            DataOutputs(fantom::DataOutputs::Control &control)
                    : DataAlgorithm::DataOutputs(control) {
                add< const Grid< 3 > >( "grid" );
                addBundle< FunctionBase >( "ColorGrid" );
                //addGraphics("streamLine");
            }
        };


            Masterarbeit( InitData& data )
             : DataAlgorithm( data ){
            }

            virtual void execute(const Algorithm::Options &options, const volatile bool &abortFlag) override{
                debugLog()<<"Der Algorithmus startet \n";

                //Setting input field and function
                std::shared_ptr<const Field<3, Matrix3> > field = options.get<Field<3, Matrix3> >("Field");
                std::shared_ptr<const Function<Matrix3> > function = options.get<Function<Matrix3> >("Field");

                // if there is no input, do nothing
                if (!field) {
                    debugLog() << "Input Field not set." << std::endl;
                    return;
                }
                // sanity check that interpolated fields really use the correct grid type. This should never fail
                std::shared_ptr<const Grid<3> > grid = std::dynamic_pointer_cast<const Grid<3> >(function->domain());
                    if (!grid) {
                    throw std::logic_error("Wrong type of grid!");
                }

                auto evaluator = field->makeEvaluator();

                std::shared_ptr<const LineSet<3> > lineSet = options.get<LineSet <3>>("LineSet");
                size_t numLines = lineSet->numLines();

                double scale= options.get<double>("Scale");
                double minValue = options.get<double>("MinValue");
                debugLog()<<"Feld und so ist gesetzt";
                //The Vector for the Points
                std::vector<Point <3>> v;
                //The Vector for the Values on the Points
                std::vector<Tensor<double>>values;

                //The Vector for the indeces
                std::vector<size_t> indices;

                //A counter for the vertices
                int verticeCounter=0;

                 debugLog() << "Beginne Schleife"<< std::endl;
                //Die Rechtecke um die Linien bauen
                for(int i=0;i<numLines;i++){
                    //debugLog() << "Durchlauf startet"<< std::endl;
                    //The indices for the Points on line i
                    const std::vector< size_t > pointIndices=lineSet->getLine(i);
                    //Calculate the first point of the Line manually to have a starting point in the Loop
                    Point3 start = lineSet->getPointOnLine(i, 0);
                    //Sanity Test
                    if(evaluator->reset(start, 0)){
                        Matrix3 startMatrix(evaluator->value());
                        auto eigensystemStart = math::getEigensystemSymmetric( startMatrix );

                        double eigenVec1 = (eigensystemStart.first[1]>minValue?eigensystemStart.first[1]:minValue);
                        double eigenVec0 = (eigensystemStart.first[0]>minValue?eigensystemStart.first[0]:minValue);

                        auto medianVecStart = scaleVector(eigensystemStart.second[1], eigenVec1*scale);
                        auto minorVecStart = scaleVector(eigensystemStart.second[0], eigenVec0*scale);
                        //auto medianVecStart = eigensystemStart.second[1];
                        //auto minorVecStart = eigensystemStart.second[0];
                        v.push_back(start+minorVecStart);
                        v.push_back(start+medianVecStart);
                        v.push_back(start-minorVecStart);
                        v.push_back(start-medianVecStart);
                        values.push_back(Tensor<double>(eigensystemStart.first[2]));
                        values.push_back(Tensor<double>(eigensystemStart.first[2]));
                        values.push_back(Tensor<double>(eigensystemStart.first[2]));
                        values.push_back(Tensor<double>(eigensystemStart.first[2]));
                        //The cap on the line
                        indices.push_back(verticeCounter);
                        indices.push_back(verticeCounter+1);
                        indices.push_back(verticeCounter+2);
                        indices.push_back(verticeCounter+3);

                        verticeCounter=verticeCounter+4;
                    } else{
                          debugLog() << "Linie konnte nicht gezeichnet werden"<< std::endl;
                    }
                    for(int j=1;j<pointIndices.size();j++){
                        //debugLog() << "Inner Loop starts"<< std::endl;
                        Point3 end = lineSet->getPointOnLine(i,j);//lineSet->getPoint(pointIndices[j]);
                        if(evaluator->reset(end, 0)){
                            Matrix3 endMatrix(evaluator->value());
                            auto eigensystemEnd = math::getEigensystemSymmetric( endMatrix );
                            double eigenVec1 = (eigensystemEnd.first[1]>minValue?eigensystemEnd.first[1]:minValue);
                            double eigenVec0 = (eigensystemEnd.first[0]>minValue?eigensystemEnd.first[0]:minValue);

                            auto medianVecEnd = scaleVector(eigensystemEnd.second[1], eigenVec1*scale);
                            auto minorVecEnd = scaleVector(eigensystemEnd.second[0], eigenVec0*scale);
                            //auto medianVecEnd = eigensystemEnd.second[1];
                            //auto minorVecEnd = eigensystemEnd.second[0];
                            
                            //verticeCounter to VerticeCounter+3 the Values of this iteration
                            v.push_back(end+minorVecEnd);
                            v.push_back(end+medianVecEnd);
                            v.push_back(end-minorVecEnd);
                            v.push_back(end-medianVecEnd);
                            values.push_back(Tensor<double>(eigensystemEnd.first[2]));
                            values.push_back(Tensor<double>(eigensystemEnd.first[2]));
                            values.push_back(Tensor<double>(eigensystemEnd.first[2]));
                            values.push_back(Tensor<double>(eigensystemEnd.first[2]));

                            //First quad
                            indices.push_back(verticeCounter);
                            indices.push_back(verticeCounter+1);
                            indices.push_back(verticeCounter-3);
                            indices.push_back(verticeCounter-4);
                            //Second quad
                            indices.push_back(verticeCounter+1);
                            indices.push_back(verticeCounter+2);
                            indices.push_back(verticeCounter-2);
                            indices.push_back(verticeCounter-3);
                            //Third quad
                            indices.push_back(verticeCounter+2);
                            indices.push_back(verticeCounter+3);
                            indices.push_back(verticeCounter-1);
                            indices.push_back(verticeCounter-2);
                            //Fourth quad
                            indices.push_back(verticeCounter+3);
                            indices.push_back(verticeCounter);
                            indices.push_back(verticeCounter-4);
                            indices.push_back(verticeCounter-1);

                            verticeCounter=verticeCounter+4;
                        }

                        //The cap on the line
                        indices.push_back(verticeCounter-1);
                        indices.push_back(verticeCounter-2);
                        indices.push_back(verticeCounter-3);
                        indices.push_back(verticeCounter-4);
                    }
                }

                const std::pair<Cell::Type, size_t> countArray[1]={
                    std::make_pair(Cell::Type::QUAD, indices.size()/4)
                };

                debugLog()<<"Schleife ist durch \n";
                std::shared_ptr< const Grid< 3 > > outputGrid = DomainFactory::makeGrid(v,(size_t) 1,countArray,indices);
                setResult("grid",outputGrid);
                debugLog()<<"Grid ist erstellt \n";
                auto outputBundle = std::make_shared< DataObjectBundle >();
                outputBundle->addContent(
                addData( outputGrid, PointSetBase::Points, std::move( values ), Precision::FLOAT64, {0} ) );//FLOAT64
                setResult( "ColorGrid", outputBundle ); 


                
               
            }

            Vector<3> scaleVector(Vector<3> v, double scale) {
                return Vector<3>(v[0] * scale, v[1] * scale, v[2] * scale);
            }

            /**
             * Converts the Quad-Indices to Triangle Indices
             * @param quad The Quad Indice Vector
             * @return the Triangle Indice Vector
             */
            std::vector<unsigned int> quadToTriangle(std::vector<size_t> quad){
                std::vector<unsigned int> triangle;
                for(int i=0;i<quad.size();i=i+4){
                    triangle.push_back(quad[i]);
                    triangle.push_back(quad[i+1]);
                    triangle.push_back(quad[i+2]);

                    triangle.push_back(quad[i+1]);
                    triangle.push_back(quad[i+2]);
                    triangle.push_back(quad[i+3]);
                }

                return triangle;
            }

            std::vector< VectorF< 3 > > convertToFloat(std::vector<Point <3>> old){
                std::vector< VectorF< 3 > > neu;
                for(int i=0;i<old.size();i++){
                    neu.push_back(VectorF<3>((float) old[i][0], (float) old[i][1], (float) old[i][2]));
                }
                return neu;
            }

            /**
             * Gets the Points on the outline of the Superquadric shape
             * @param numPoints the Number of Points
             * @param eigenSystem the eigenSystem, i.e. the eigenvalues and eigenvectors
             * @return a vector of the Points on the outer edge of the Superquadric shape*
             */
            std::vector<Vector<3>> getPoints(int numPoints, std::pair< std::array< double, 3 >, std::array< Tensor< double, 3 >, 3 > > eigenSystem){
                
                double sum = eigenSystem.first[0]+eigenSystem.first[1]+eigenSystem.first[2];
                double cl=(eigenSystem.first[2]-eigenSystem.first[1])/sum;
                double cp=(eigenSystem.first[1]-eigenSystem.first[0])/sum;

                double gamma = 1.0;//Am Ende vom Nutzer eingeben lassen
                double alpha;
                if(cl>=cp){
                    alpha=pow(1-cp,gamma);
                }else{
                    alpha=pow(1-cl,gamma);
                }

                double pifraction=M_PI/numPoints;
                std::vector<Vector<3>> points;
                //Initializing the Points. The center is (0,0,0)
                for(int i=0;i<=numPoints;i++){
                    double costheta=cos(pifraction*i);
                    double sintheta=sin(pifraction*i);
                    double x = sign(costheta)*pow(abs(costheta),alpha);
                    double y = sign(sintheta)*pow(abs(sintheta),alpha);
                    points.push_back(Vector<3>(x,y,0.0));
                }
                return points;
            }

            double sign(double a){
                if(a<0){
                    return -1.0;
                } else{
                    return 1.0;
                }
            }

            //auto eigensystem = math::getEigensystemSymmetric( tensor\ );
            /*auto majorVec = eigensystem.second[2];
            auto medianVec = eigensystem.second[1];
            auto minorVec = eigensystem.second[0];
            auto majorVal = eigensystem.first[2];
            auto medianVal = eigensystem.first[1];
            auto minorVal = eigensystem.first[0];*/

            // auto lineset = std::make_shared< LineSet >();
            // std::vector<size_t> line;
            //line.push_back( lineset->addPoint( points[i] ) );
            //lineset->addLine( line );
            //line.clear();
 
    };
    
    AlgorithmRegister< Masterarbeit > dummy( "Masterarbeit/Visualisierung", "Visualisierung" );
}
