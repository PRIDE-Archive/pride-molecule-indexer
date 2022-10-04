package uk.ac.ebi.pride.archive.indexer.utility;

import org.ehcache.Cache;
import org.ehcache.CacheManager;
import org.ehcache.config.Configuration;
import org.ehcache.config.builders.CacheConfigurationBuilder;
import org.ehcache.config.builders.CacheManagerBuilder;
import org.ehcache.config.builders.ResourcePoolsBuilder;
import org.ehcache.xml.XmlConfiguration;

import java.net.URL;
import java.util.List;

public class AppCacheManager {
    private static AppCacheManager instance = null;
    private CacheManager cacheManage = null;

    private static String INPUT_SPECTRA_CACHE = "InputSpectra";
    private static String CLUSTERS_CACHE = "Clusters";
    private static String PROTEIN_TO_PSMS_CACHE = "ProteinToPsms";

    private AppCacheManager() {

        cacheManage = CacheManagerBuilder.newCacheManagerBuilder()
                .withCache(INPUT_SPECTRA_CACHE, CacheConfigurationBuilder
                        .newCacheConfigurationBuilder(String.class, Long.class, ResourcePoolsBuilder.heap(100000)))
                .withCache(CLUSTERS_CACHE, CacheConfigurationBuilder
                        .newCacheConfigurationBuilder(Integer.class, Integer.class, ResourcePoolsBuilder.heap(10000)))
                .withCache(PROTEIN_TO_PSMS_CACHE, CacheConfigurationBuilder
                        .newCacheConfigurationBuilder(String.class, List.class, ResourcePoolsBuilder.heap(10000)))
                .build();
        cacheManage.init();

        System.out.println("Cache Initialized");
    }

    public static AppCacheManager getInstance() {
        if (instance == null) {
            instance = new AppCacheManager();
        }
        return instance;
    }

    public Cache<Integer, Integer> getClustersCache(){
        return this.cacheManage.getCache(CLUSTERS_CACHE, Integer.class, Integer.class);
    }

    public Cache<String, Long> getPrideJsonSpectra() {
        return this.cacheManage.getCache(INPUT_SPECTRA_CACHE, String.class, Long.class);
    }


    public Cache<String, ? extends List> getProteinToPsmsCache() {
        return this.cacheManage.getCache(PROTEIN_TO_PSMS_CACHE, String.class, List.class);
    }
}
