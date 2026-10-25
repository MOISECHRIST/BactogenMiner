#!/usr/bin/env nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
UNTAR_DB : Untar databases archives when they are given as .tar.gz files 
    INPUT:
        archive : Path to the database archive (.tar.gz)
    OUTPUT :
        db : Path to the database stored in params.db_cache
*/

process UNTAR_DB {
    storeDir params.db_cache

    input:
        path archive

    output:
        path "${archive.name.replaceAll(/\.tar\.gz.*$/, '')}", emit: db

    script:
    def name = archive.name.replaceAll(/\.tar\.gz.*$/, '')
    """
    mkdir ${name}
    tar -xzf ${archive} -C ${name} --strip-components=1
    """
}


def isArchive(path) {
    def string = path.toString()
    return string.endsWith('.tar.gz') || string.contains('.tar.gz?')
}

workflow PREPARE_DB {
    take:
        db_path   // channel: value, local folder or path/URL to a .tar.gz

    main:
        routed = db_path.branch { p ->
            archive: isArchive(p)
            folder:  true// single item, reused by every sample
        }

        UNTAR_DB(routed.archive.map { p -> file(p) })

        db_ch = UNTAR_DB.out.db
            .mix(routed.folder.map { p -> file(p) })
            .first()

    emit:
        db = db_ch
}

process GET_DB {
    tag "${name}"
    storeDir "${params.db_cache}/archives"

    input:
        tuple val(name), val(url), val(md5)

    output:
        path "${name}.tar.gz", emit: db

    script:
    def expected = (md5 ?: '').toString().replaceFirst(/^md5:/, '')
    """
    wget --quiet --tries=3 --timeout=60 -O ${name}.tar.gz "${url}"

    if [ -n "${expected}" ]; then
        echo "${expected}  ${name}.tar.gz" | md5sum -c -
    fi

    tar -tzf ${name}.tar.gz > /dev/null
    """
}