void MypolysConstructor(int rank,
                        const std::vector<Polygon*>& polygons,
                        std::vector<int>& mypolys)
{
    for (size_t i = 0; i < polygons.size(); ++i) 
        if (polygons[i]->getRank() == rank) {
            mypolys.push_back(i);
        }

    for (size_t i = 0; i < facets.size(); ++i) 
        for (size_t j = 0; j < facets[i]->getNeigbors().size(); ++j) 
            if (cells[facets[i]->getNeigbors()[j]]->getRank() == rank) {
                myfacets.push_back(i);
                break;
            }
}

void UnstructMesh::setMpiRank(int rank) 
{
    MypolysConstructor(rank, cells, mycells);
}

