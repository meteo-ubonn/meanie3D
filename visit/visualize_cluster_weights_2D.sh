
if [ "X$1" = "X" ]
then
echo "visualize_cluster_weights_3D.sh <source file> <variable>"
echo "Creates a python script for cluster creation and tracking and runs it in Visit"
exit 0
fi

if [ "X$2" = "X" ]
then
echo "visualize_cluster_weights_3D.sh <source file> <variable>"
echo "Creates a python script for cluster creation and tracking and runs it in Visit"
exit 0

fi

if [ "X${VISIT_EXECUTABLE}" = "X" ]
then
echo "Please set environment variable VISIT_EXECUTABLE"
exit 0
fi
if [ "X${MEANIE3D_HOME}" = "X" ]
then
echo "Please set environment variable MEANIE3D_HOME"
exit 0
fi

#DL_PATH=$MEANIE3D_HOME/Release
DL_PATH=/usr/local/lib

SCRIPTFILE="/tmp/visit-$RANDOM.py"
ESCAPED_SOURCE_DIR=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")
ESCAPED_DL_PATH=$(echo $DL_PATH | sed -e "s/\//\\\\\//g")

cat $MEANIE3D_HOME/visit/visualize_cluster_weights_2D.py | sed -e "s/SOURCE_FILE_P/$ESCAPED_SOURCE_DIR/g" | sed -e "s/M3D_HOME_P/$ESCAPED_MEANIE3D_HOME/g" | sed -e "s/VAR_NAME_P/$2/g"> $SCRIPTFILE
#${VISIT_EXECUTABLE} -cli -nowin -s $SCRIPTFILE
${VISIT_EXECUTABLE} -s $SCRIPTFILE